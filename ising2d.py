import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from scipy.fftpack import fft, ifft, ifftshift
from scipy.optimize import curve_fit



class ising2d():
    def __init__(self, L, T, B, algorithm='metropolis'):
        self.T = T
        self.B = B
        self.L = L
        self.N = L**2
        self.state = np.random.choice([-1,1], size=(L,L))
        
        self.algorithm = algorithm
        self.equilibrium = False

        
        self.corrtime = None
        self.energy_evolution = None
        self.autocorrelation = None
        self.delays = None
        self.observables = []

        self.__energy()
        self.__magnetization()
        self.__probability()


    def update_system(self, T, B):
        self.T = T
        self.B = B
        self.equilibrium = False
        self.correlation_time = None
        self.energy_evolution = None
        self.autocorrelation = None
        self.delays = None
        self.observables = []

        self.E = self.__energy()
        self.M = self.__magnetization()
        self.__probability()
    
    def update_microstate(self):
        """ Flip spins until the energy correlations are gone and an independent configuration is generated """
        if self.equilibrium:
            if self.corrtime:
                    self.__spinflip(5*self.corrtime)
                    self.__save_observables()
            else:
                raise RuntimeError('The correlation time has not been set') 
        else:
            raise RuntimeError('The system is not in equilibrium')

    
    def thermalize(self):
        """ Perform enough spin flip operations that the system reaches thermal equilibrium """
        if self.algorithm == 'metropolis':
            steps = 100*self.N**2
        else:
            steps = 100*self.N
        self.__spinflip(steps)
        self.equilibrium = True

    def correlation_time(self, plot=False):
        """ Flip spins and keep track of energy evolution over time to collect correlation data """
        self.__energy_evolution()
        self.__autocorrelation()
        self.delays = np.arange(len(self.autocorrelation))
        popt, pcov = curve_fit(self.__exponential, self.delays, self.autocorrelation, p0=[self.N])
        self.corrtime = popt[0]
        if plot:
            pl.plot(self.delays, self.autocorrelation, label='Autocorrelation of Energy')
            pl.plot(self.delays, self.__exponential(self.delays, self.corrtime), label='Single Exponential Fit')
            pl.legend(loc='best')
            pl.show()

    ##private internal utility functions
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
            neighbours = np.array([self.state[i, (j+1)%self.L], self.state[i, (j-1)%self.L], self.state[(i+1)%self.L, j], self.state[(i-1)%self.L, j]],dtype=np.int64)
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
        dE = 0
        dM = 0
        pass


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
            self.__spinflip(100*self.N**2, save=True)
        elif self.algorithm == 'wolff':
            raise NotImplementedError('wolff not yet implemented')

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
        row = {'L': self.L, 'N': self.N, 'T': self.T, 'B': self.B, 'E': self.E, 'M': self.M, 'correlation_time': self.corrtime}
        self.observables.append(row)
        
    ## output functions
    def print_energy_evolution(self, filename):
        """ Save the time evolution of the energy to a csv file """
        np.savetxt(filename, self.energy_evolution, delimiter=',')
    
    def print_state(self, filename):
        """ Print a 2D binary matrix representing the spins in the system """
        np.savetxt(filename, self.state, delimiter=',')

    def print_observables(self, filename):
        """ Save all of the generated observables in a csv file """
        pd.DataFrame(self.observables).to_csv(filename, sep=',', index=False)

    def print_autocorrelation(self, filename):
        """ Save the autocorrelation function in a csv file """
        np.savetxt(filename, self.autocorrelation, delimiter=',')


