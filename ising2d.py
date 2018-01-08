import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from scipy.fftpack import fft, ifft, ifftshift
from scipy.optimize import curve_fit
from tqdm import tqdm
import itertools
from collections import deque
import os

class ising2d():
    def __init__(self, temperatures, fields, sizes, microstates, algorithm='metropolis', output_folder='.', save_states = 0):
        self.algorithm = algorithm
        self.observables = []
        self.correlations = []
        self.output_folder = output_folder
        self.temperatures = temperatures
        self.fields = fields
        self.sizes = sizes
        if algorithm == 'wolff':
            self.eqbm_window = 50
            self.eqbm_zerocount = 20
        elif algorithm == 'metropolis':
            self.eqbm_window = 500
            self.eqbm_zerocount = 200
        self.microstates = microstates
        self.save_states = save_states
        self.saved_states = 0
        

        if any(np.array(temperatures) < 1):
            raise ValueError('The Monte Carlo method cannot be reliably used for T < 1')
        if any(np.absolute(fields) > 0.1) and algorithm == 'wolff':
            raise ValueError('The Wolff Algorithm performs poorly for B > 0.1, consider using a smaller field or use Metropolis if a large field is needed')

        paths = [output_folder, output_folder+'/states', output_folder+'/correlations']
        for path in paths:
            try: 
                os.makedirs(path)
            except OSError:
                if not os.path.isdir(path):
                    raise

    def run(self):
        """ generate mictostates for all possible ensembles generated by the lists of size, temperature, and magnetic field """
        ensembles = itertools.product(self.sizes, self.temperatures, self.fields)
        pbar = tqdm(ensembles, total=len(self.sizes)*len(self.temperatures)*len(self.fields))
        for ensemble in pbar:
            L = ensemble[0]
            T = ensemble[1]
            B = ensemble[2]
            pbar.set_description('(L,T,B) = ({0}, {1:.3g}, {2:.3g})'.format(L,T,B))
            self._update_system(L,T,B)
            self.correlation_length = np.zeros(self.L/2)
            self._print_energy_evolution()
            self._print_autocorrelation()
            for k in tqdm(range(self.microstates), desc = 'Production'):
                self._update_microstate()
                if self.saved_states < self.save_states:
                    self._save_state()
                self._correlation_length()
            self.correlation_length /= float(self.microstates)
            self._fit_correlation_length()
            self._print_correlation_length()
            self._save_correlations()
        self._print_observables()
        self._print_correlations()

    ##private internal utility functions
    
    def _thermalize(self):
        """ Perform enough spin flip operations that the system reaches thermal equilibrium """
        self.zerocount = 0
        self.thermalsteps = 0
        if self.algorithm == 'metropolis':
            steps = self.N**2
        else:
            steps = np.maximum(self.N/10, 500)
        while self.zerocount < self.eqbm_zerocount:
            self._spinflip(steps, mode='Thermalize')
            self.thermalsteps += steps
        
    def _update_system(self, L, T, B):
        """ set a new ensemble and equilibrate it """
        self.T = T
        self.B = B
        self.L = L
        self.N = L**2
        self.state = np.random.choice([-1,1], size=(L,L))
        self.state = self.state.astype(np.int64)
        self.corrtime = None
        self.energy_evolution = None
        self.autocorrelation = None
        self.delays = None
        self.saved_states = 0
        self._energy()
        self._magnetization()
        self._probability()
        self._thermalize()
        self._correlation_time()
    
    def _update_microstate(self):
        """ Flip spins until the energy correlations are gone and an independent configuration is generated """
        self._spinflip(5*int(self.corrtime+1), mode = 'Production')
        self._save_observables()

    def _correlation_time(self):
        """ Flip spins and keep track of energy evolution over time to collect correlation data """
        self._energy_evolution()
        self._autocorrelation()
        self.delays = np.arange(len(self.autocorrelation))
        if self.algorithm == 'metropolis':
            default = self.N
        else:
            default = 3
        p0 = [next((i for i in self.delays if self.autocorrelation[i] < 0), default)/3.0]
        popt, pcov = curve_fit(self._exponential, self.delays, self.autocorrelation, p0=p0)
        self.corrtime = popt[0]
            
    def _spinflip(self, steps, mode=None):
        """ perform a single spin update step using the given algorithm """
        if self.algorithm == 'metropolis':
            self._metropolis(steps, mode)
        elif self.algorithm == 'wolff':
            self._wolff(steps, mode)
        else:
            raise NotImplementedError('The {0} algorithm is not supported'.format(self.algorithm))

    def _metropolis(self, steps, mode):
        """ perform spin update steps using the Metropolis algorithm """

        label = mode
        if mode == 'Thermalize':
            save = False
        elif mode == 'Autocorrelation':
            save = True
        elif mode == 'Production':
            save = False
            label = None
        else:
            raise NotImplementedError('{0} is not a valid mode'.format(mode))
        
        if mode=='Thermalize':
            self.energy_sign = deque()
        if save:
            self.energy_evolution = np.zeros(steps)
        spins = np.random.randint(0, self.L, size=(steps, 2))
        if label is not None:
            iterator = tqdm(spins, desc = label)
        else:
            iterator = spins
        step = 0
        for spin in iterator:
            i = spin[0] % self.L
            j = spin[1] % self.L
            s = self.state[i,j]
            ss = (s+1)//2
            neighbours = np.array([self.state[i, (j+1)%self.L], self.state[i, j-1], self.state[(i+1)%self.L, j], self.state[i-1, j]],dtype=np.int64)
            upneighbours = np.sum((neighbours+1)//2)
            p = self.probability[ss, upneighbours]
            flipped = False
            if p >= 1 or np.random.rand() < p:
                self.state[i,j] *= -1
                self.E += self.energytable[ss, upneighbours]
                self.M += 2*self.state[i,j]
                flipped = True
            if save:
                self.energy_evolution[step] = self.E
                step += 1
            if mode=='Thermalize':
                added = False
                if flipped:
                    self.energy_sign.append(np.sign(self.energytable[ss, upneighbours]))
                    added = True
                if len(self.energy_sign) > self.eqbm_window:
                    self.energy_sign.popleft()
                if np.sum(self.energy_sign) == 0 and len(self.energy_sign) == self.eqbm_window and added:
                    self.zerocount += 1
                if self.zerocount == self.eqbm_zerocount:
                    break

    def _wolff(self, steps, mode):
        """ perform a spin cluster update step using the Wolff algorithm """
        label = mode
        if mode == 'Thermalize':
            save = False
        elif mode == 'Autocorrelation':
            save = True
        elif mode == 'Production':
            save = False
            label=None
        else:
            raise NotImplementedError('{0} is not a valid mode'.format(mode))
        

        
        if mode=='Thermalize':
            self.energy_sign = deque()
        if save:
            self.energy_evolution = np.zeros(steps)
            
        if label:
            iterator = tqdm(range(steps), desc=label)
        else:
            iterator = range(steps)
        for k in iterator:
            cluster, sign = self._build_cluster(self.probability)
            oldE = self.E
            if self.B == 0 or (np.sign(self.B) == -sign) or (np.random.rand() < np.exp(-2*sign*self.B*np.sum(cluster))):
                self.state[cluster == 1] *= -1
                self._energy()
                self.M -= 2*np.sum(cluster)*sign
            dE = oldE - self.E
            if save:
                self.energy_evolution[k] = self.E
            if mode=='Thermalize':
                added = False
                if dE != 0:
                    self.energy_sign.append(np.sign(dE))
                    added = True
                if len(self.energy_sign) > self.eqbm_window:
                    self.energy_sign.popleft()
                if np.sum(self.energy_sign) == 0 and len(self.energy_sign) == self.eqbm_window and added:
                    self.zerocount += 1
                if self.zerocount == self.eqbm_zerocount:
                    break
                    
            

    def _build_cluster(self, prob, seed=None):
        """ build a cluster of like-spin nearest neighbours to be flipped all at once"""
        cluster = np.zeros(self.state.shape, dtype=np.int64)
        pocket = []
        if seed == None:
            seed = np.squeeze(np.random.randint(0, self.L, size=(1,2)).tolist())
        sign = self.state[seed[0], seed[1]]
        cluster[seed[0], seed[1]] = 1
        pocket.append(seed)
        pocketnum = 1
        index = 0
        while index < pocketnum:
            i = pocket[index][0]
            j = pocket[index][1]
            neighbours = [[(i+1)%self.L, j], [(i-1)%self.L, j], [i,(j+1)%self.L], [i, (j-1)%self.L]]
            for neighbour in neighbours:
                x = neighbour[0]
                y = neighbour[1]
                if self.state[i,j] == self.state[x,y] and cluster[x,y] != 1 and np.random.rand() < prob:
                    pocket.append([x,y])
                    cluster[x,y] = 1
                    pocketnum += 1
            index += 1
        return cluster, sign
                    
    def _exponential(self, n, n0):
        """ return a single exponential function for fitting """
        return np.exp(-(n/n0))

    def _offset_exponential(self, n, n0, a, p):
        """ return a single exponential function for fitting """
        return a*np.exp(-(n/n0))/n**p

    def _energy(self):
        """ Calculate the total energy of the system """
        self.E = -np.sum(self.state*(np.roll(self.state, 1, axis=0) + np.roll(self.state, 1, axis=1) + self.B))
            

    def _magnetization(self):
        """ Calculate the total magnetization of the system """
        self.M = np.sum(self.state)

    def _autocorrelation(self):
        """ Calculate the autocorrelation of the energy of the system using that fact that the autocorrelation is the Fourier Transform of the PSD """
        energy = self.energy_evolution
        maxdelay = len(energy)/5
        xp = ifftshift((energy - np.average(energy))/np.std(energy))
        n = len(xp)
        xp = np.r_[xp[:n//2], np.zeros_like(xp), xp[n//2:]]
        f = fft(xp)
        S = np.absolute(f)**2
        R = ifft(S)
        self.autocorrelation = (np.real(R)[:n//2]/(np.arange(n//2)[::-1]+n//2))[:maxdelay]

    def _correlation_length(self):
        pass
##        maxlength = self.L/2
##        mean = np.average(self.state)
##        std = np.std(self.state)
##        correlation = np.zeros(maxlength, dtype=np.float64)
##        samples = np.zeros(maxlength, dtype=np.float64)
##        for i in range(self.L):
##            for j in range(self.L):
##                for k in range(maxlength):
##                    correlation[k] += (self.state[i,j]-mean)/std * ((self.state[i, (j-k)%self.L]-mean) + (self.state[i, (j+k)%self.L]-mean) + (self.state[(i+k)%self.L, j]-mean) + (self.state[(i-k)%self.L, j]-mean))/std
##                    samples[k] += 4
##        correlation /= samples
##        self.correlation_length += correlation

    def _fit_correlation_length(self):
        self.corrlength = 0
        self.eta = 0
        
##        p0 = [5, 1, 0.25]
##        x = np.arange(len(self.correlation_length))
##        if self.T < 2/np.log(1+np.sqrt(2)):
##            self.corrlength = 0
##            self.eta = 0
##        else:
##            popt, pcov = curve_fit(self._offset_exponential, x[1:], self.correlation_length[1:], p0=p0)
##            self.corrlength = popt[0]
##            self.eta = popt[2]

    def _energy_evolution(self):
        """ Flip spins and keep track of energy evolution over time to collect correlation data """
        self._spinflip(np.maximum(1000, 2*self.thermalsteps), mode='Autocorrelation')

    def _probability(self):
        """ pre-define the spin-flip/cluster addition probabilities """
        if self.algorithm == 'metropolis':
            self.energytable = np.zeros((2,5))
            for j in range(5):
                self.energytable[0,j] = 2*(4.0 - 2.0*j - 2*self.B) #spin down, with j neighbours spin up
                self.energytable[1,j] = 2*(2.0*j - 4.0 + 2*self.B) #spin up, with j neighnours spin up
            self.probability = np.minimum(1.0, np.exp(-self.energytable/self.T))     
        elif self.algorithm == 'wolff':
            self.probability = 1.0 - np.exp(-2.0/self.T)
        else:
            raise NotImplementedError('The {0} algorithm is not supported'.format(self.algorithm))
            
    def _save_observables(self):
        """ Add a row of observables to the list of saved microstates """
        row = {'L': self.L, 'N': self.N, 'T': self.T, 'B': self.B, 'E': self.E, 'M': self.M}
        self.observables.append(row)

    def _save_correlations(self):
        row = {'L': self.L, 'N': self.N, 'T': self.T, 'B': self.B, 'correlation_time': self.corrtime, 'correlation_length': self.corrlength, 'eta': self.eta}
        self.correlations.append(row)

    def _save_state(self):
        """ Save the time evolution of the energy to a csv file """
        np.savetxt(self.output_folder + '/states/state_T={0:.6g}_B={1:.6g}_L={2}_{3}.csv'.format(self.T,self.B,self.L, self.saved_states), self.state, delimiter=',')
        pl.imshow(self.state, interpolation='none')
        pl.xlabel('X Spin Index')
        pl.ylabel('Y Spin Index')
        pl.title('T={0:.6g} B={1:.6g} L={2}'.format(self.T,self.B,self.L))
        pl.savefig(self.output_folder + '/states/state_T={0:.6g}_B={1:.6g}_L={2}_{3}.png'.format(self.T,self.B,self.L, self.saved_states))
        pl.close('all')
        self.saved_states += 1
       
    def _print_energy_evolution(self):
        """ Save the time evolution of the energy to a csv file """
        np.savetxt(self.output_folder + '/correlations/energy_evolution_T={0:.6g}_B={1:.6g}_L={2}.csv'.format(self.T,self.B,self.L), self.energy_evolution, delimiter=',')

    def _print_observables(self):
        """ Save all of the generated observables in a csv file """
        pd.DataFrame(self.observables).to_csv(self.output_folder + '/observables.csv', sep=',', index=False)

    def _print_correlations(self):
        """ Save all of the generated observables in a csv file """
        pd.DataFrame(self.correlations).to_csv(self.output_folder + '/correlations.csv', sep=',', index=False)

    def _print_autocorrelation(self):
        """ Save the autocorrelation function in a csv file """
        np.savetxt(self.output_folder + '/correlations/autocorrelation_T={0:.6g}_B={1:.6g}_L={2}.csv'.format(self.T,self.B,self.L), self.autocorrelation, delimiter=',')

    def _print_correlation_length(self):
        """ Save the autocorrelation function in a csv file """
        np.savetxt(self.output_folder + '/correlations/correlation_length_T={0:.6g}_B={1:.6g}_L{2}.csv'.format(self.T,self.B,self.L), self.correlation_length, delimiter=',')


