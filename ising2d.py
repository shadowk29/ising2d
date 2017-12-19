import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from scipy.fftpack import fft, ifft, ifftshift
from scipy.optimize import curve_fit
import progressbar

class ising2d():
    def __init__(self, temperatures, fields, sizes, microstates, algorithm='metropolis', output_folder='.', verbose = False):
        self.algorithm = algorithm
        self.observables = []
        self.output_folder = output_folder
        self.temperatures = temperatures
        self.fields = fields
        self.sizes = sizes
        self.microstates = microstates
        self.verbose = verbose

        if any(temperatures < 1):
            raise ValueError('The Monte Carlo method cannot be reliably used for T < 1')
        if any(np.absolute(fields) > 0) and algorithm == 'wolff':
            raise ValueError('The Wolff Algorithm can only be used when B = 0')

    def run(self):
        done = 0
        with progressbar.ProgressBar(max_value=self.microstates * len(self.sizes) * len(self.fields) * len(self.temperatures)) as bar:
            for L in self.sizes:
                for T in self.temperatures:
                    for B in self.fields:
                        self.__update_system(L,T,B)
                        if self.verbose:
                            self.__print_energy_evolution()
                            self.__print_autocorrelation()
                        for k in range(self.microstates):
                            self.__update_microstate()
                            done += 1
                            bar.update(done)
        self.__print_observables()

    ##private internal utility functions
    
    def __thermalize(self):
        """ Perform enough spin flip operations that the system reaches thermal equilibrium """
        if self.algorithm == 'metropolis':
            steps = self.N**2
        else:
            steps = self.N
        self.__spinflip(steps)
        
    def __update_system(self, L, T, B):
        self.T = T
        self.B = B
        self.L = L
        self.N = L**2
        self.state = np.random.choice([-1,1], size=(L,L))
        self.corrtime = None
        self.energy_evolution = None
        self.autocorrelation = None
        self.delays = None
        self.__energy()
        self.__magnetization()
        self.__probability()
        self.__thermalize()
        self.__correlation_time()
    
    def __update_microstate(self):
        """ Flip spins until the energy correlations are gone and an independent configuration is generated """
        if self.ready:
            self.__spinflip(5*self.corrtime)
            self.__save_observables()
        else:
            raise RuntimeError('The ensemble (L,B,T) has not been specified')

    def __correlation_time(self, plot=False):
        """ Flip spins and keep track of energy evolution over time to collect correlation data """
        self.__energy_evolution()
        self.__autocorrelation()
        self.delays = np.arange(len(self.autocorrelation))
        popt, pcov = curve_fit(self.__exponential, self.delays, self.autocorrelation, p0=[self.N])
        self.corrtime = int(popt[0])
        if self.corrtime == 0:
            self.corrtime = 1
        if plot:
            pl.plot(self.delays, self.autocorrelation, label='Autocorrelation of Energy')
            pl.plot(self.delays, self.__exponential(self.delays, self.corrtime), label='Single Exponential Fit')
            pl.legend(loc='best')
            pl.show()
            
    def __spinflip(self, steps, save=False):
        """ perform a single spin update step using the given algorithm """
        if self.algorithm == 'metropolis':
            self.__metropolis(steps, save)
        elif self.algorithm == 'wolff':
            self.__wolff(steps, save)
        else:
            raise NotImplementedError('The {0} algorithm is not supported'.format(self.algorithm))

    def __metropolis(self, steps, save):
        """ perform spin update steps using the Metropolis algorithm """
        if save:
            self.energy_evolution = np.zeros(steps)
        spins = np.random.randint(0, self.L, size=(steps, 2))
        step = 0
        for spin in spins:
            i = spin[0] % self.L
            j = spin[1] % self.L
            s = self.state[i,j]
            ss = (s+1)//2
            neighbours = np.array([self.state[i, (j+1)%self.L], self.state[i, j-1], self.state[(i+1)%self.L, j], self.state[i-1, j]],dtype=np.int64)
            upneighbours = np.sum((neighbours+1)//2)
            p = self.probability[ss, upneighbours]
            if p >= 1 or np.random.rand() < p:
                self.state[i,j] *= -1
                self.E += self.energytable[ss, upneighbours]
                self.M += 2*self.state[i,j]
            if save:
                self.energy_evolution[step] = self.E
                step += 1

    def __wolff(self, steps, save):
        """ perform a spin cluster update step using the Wolff algorithm """
        if save:
            self.energy_evolution = np.zeros(steps)
        for k in range(steps):
            cluster, sign = self.__build_cluster(self.probability)
            self.state[cluster == 1] *= -1
            self.__energy()
            self.M -= 2*np.sum(cluster)*sign
            if save:
                self.energy_evolution[k] = self.E
        

    def __build_cluster(self, prob, seed=None):
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

    def __domain_size(self):
        clusters = np.zeros(self.state.shape, dtype=np.int64)
        clustersize = 0
        maxclustersize = 0
        count = 0
        for i in range(self.L):
            for j in range(self.L):
                if clusters[i,j] == 0:
                    cluster, sign = self.__build_cluster(1.0, [i,j])
                    clustersize = np.sum(cluster)
                    if clustersize > maxclustersize:
                        maxclustersize = clustersize
                    count += 1
                    clusters += cluster
        self.cluster_count = count
        self.maxclustersize = maxclustersize
        print np.sum(clusters)
                    
    def __exponential(self, n, n0):
        """ return a single exponential function for fitting """
        return np.exp(-(n/n0))

    def __energy(self):
        """ Calculate the total energy of the system """
        self.E = -np.sum(self.state*(np.roll(self.state, 1, axis=0) + np.roll(self.state, 1, axis=1) + self.B))
            

    def __magnetization(self):
        """ Calculate the total magnetization of the system """
        self.M = np.sum(self.state)

    def __autocorrelation(self):
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

    def __energy_evolution(self):
        """ Flip spins and keep track of energy evolution over time to collect correlation data """
        if self.algorithm == 'metropolis':
            self.__spinflip(50*self.N**2, save=True)
        elif self.algorithm == 'wolff':
            self.__spinflip(50*self.N, save=True)

    def __probability(self):
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
            
    def __save_observables(self):
        """ Add a row of observables to the list of saved microstates """
        self.__domain_size()
        row = {'L': self.L, 'N': self.N, 'T': self.T, 'B': self.B, 'E': self.E, 'M': self.M, 'correlation_time': self.corrtime, 'cluster_count': self.cluster_count, 'max_cluster_size': self.maxclustersize}
        self.observables.append(row)
        
    def __print_energy_evolution(self):
        """ Save the time evolution of the energy to a csv file """
        np.savetxt(self.output_folder + '/energy_evolution_T={0:.2f}_B={1:.2f}_L{2}.csv'.format(self.T,self.B,self.L), self.energy_evolution, delimiter=',')

    def __print_observables(self):
        """ Save all of the generated observables in a csv file """
        pd.DataFrame(self.observables).to_csv(self.output_folder + '/observables.csv'.format(self.T,self.B,self.L), sep=',', index=False)

    def __print_autocorrelation(self):
        """ Save the autocorrelation function in a csv file """
        np.savetxt(self.output_folder + '/autocorrelation_T={0:.2f}_B={1:.2f}_L{2}.csv'.format(self.T,self.B,self.L), self.autocorrelation, delimiter=',')


