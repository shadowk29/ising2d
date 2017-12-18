import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from scipy.fftpack import fft, ifft, ifftshift
from scipy.optimize import curve_fit



class ising2d():
    def __init__(L, T, B, algorithm='metropolis'):
        self.T = T
        self.B = B
        self.L = L
        self.N = L**2
        self.state = np.random.choice([-1,1], size=(L,L))
        self.E = self.__energy()
        self.M = self.__magnetization()
        self.algorithm = algorithm
        self.equilibrium = False

        
        self.correlation_time = None
        self.energy_evolution = None
        self.autocorrelation = None
        self.delays = None
        self.observables = []

        self.__probability()


    def update_system(self, T, B):
        self.T = T
        self.B = B
        self.equilibrium = False
        self.correlation_time = None
    
    def update_microstate(self):
        """ Flip spins until the energy correlations are gone and an independent configuration is generated """
        if self.equilibrium:
            if self.correlation_time:
                for i in range(3*correlation_time):
                    self.__spinflip()
            else:
                raise RuntimeError('The correlation time has not been set') 
        else:
            raise RuntimeError('The system is not in equilibrium')

    
    def thermalize(self):
        """ Perform enough spin flip operations that the system reaches thermal equilibrium """
        if self.algorithm == 'metropolis':
            steps = N**2
        else:
            steps = N
        for i in range(steps):
            dE, dM = self.__spinflip()
            self.E += dE
            self.M += dM
        self.equilibrium = True

    def correlation_time(self, plot=False):
        """ Flip spins and keep track of energy evolution over time to collect correlation data """
        self.delays = np.arange(len(self.autocorrelation))
        popt, pcov = curve_fit(exponential, delays, self.autocorrelation, p0=[self.N])
        self.correlation_time = popt[0]
        if plot:
            pl.plot(self.delays, self.autocorrelation, label='Autocorrelation of Energy')
            pl.plot(self.delays, self.exponential(self.delays, n0), label='Single Exponential Fit')
            pl.legend(loc='best')
            pl.show()

    ##private internal utility functions
    def __spinflip(self):
        """ perform a single spin update step using the given algorithm """
        if self.algorithm = 'metropolis':
            dE, dM = self.__metropolis()
        elif self.algorithm = 'wolff':
            dE, dM = self.__wolff()
        else:
            raise NotImplementedError('The {0} algorithm is not supported'.format(self.algorithm))
        return dE, dM

    def __metropolis(self):
        """ perform a single spin update step using the Metropolis algorithm """
        pass

    def __wolff(self):
        """ perform a spin cluster update step using the Wolff algorithm """
        pass


    def __exponential(self, n, n0):
        """ return a single exponential function for fitting """
        return np.exp(-(n/n0))

    def __energy(self):
        """ Calculate the total energy of the system """
        energy = 0
        for i in range(self.L):
            for j in range(self.L):
                energy -= self.state[i,j]*(self.state[i, (j+1)%self.L] + self.state[(i+1)%self.L, j] + self.B)
        self.E = energy

    def __magnetization(self):
        """ Calculate the total magnetization of the system """
        self.M = np.sum(state)

    def __autocorrelation(self, energy):
        """ Calculate the autocorrelation of the energy of the system using that fact that the autocorrelation is the Fourier Transform of the PSD """
        xp = ifftshift((energy - np.average(energy))/np.std(energy))
        n = len(xp)
        xp = np.r_[xp[:n/2], np.zeros_like(xp), xp[n/2:]]
        f = fft(xp)
        S = np.absolute(f)**2
        R = ifft(S)
        self.autocorrelation = np.real(R)[:n/2]/(np.arange(n/2)[::-1]+n/2)

    def __energy_evolution(self):
        """ Flip spins and keep track of energy evolution over time to collect correlation data """
        steps = 20*self.N
        self.energy_evolution = np.zeros(steps)
        for i in range(steps):
            dE, dM = self.__spinflip()
            self.E += dE
            self.M += dM
            self.energy_evolution[i] = self.E

    def __probability(self):
        if self.algorithm == 'metropolis':
            self.energytable = np.zeros((2,5))
            for j in range(5):
                self.energytable[0,j] = 4.0 - 2.0*j - 2*self.B #spin down, with j neighbours spin up
                self.energytable[1,j] = 2.0*j - 4.0 + 2*self.B #spin up, with j neighnours spin up
            self.probability = np.minimum(1.0, np.exp(-2.0/self.T*self.energytable))     
        elif self.algorithm == 'wolff':
            self.probability = 1.0 - np.exp(-2.0/self.T)
        else:
            raise NotImplementedError('The {0} algorithm is not supported'.format(self.algorithm))
            

    ## output functions
    def print_energy_evolution(self, filename):
        """ Save the time evolution of the energy to a csv file """
        pass
    
    def print_state(self, filename):
        """ Print a 2D binary matrix representing the spins in the system """
        pass

    def print_observables(self, filename):
        """ Save all of the generated observables in a csv file """
        pass

    def print_autocorrelation(self, filename):
        """ Save the autocorrelation function in a csv file """
        pass

    def save_observables(self):
        """ Add a row of observables to the list of saved microstates """
        row = {'L': self.L, 'N': self.N, 'T': self.T, 'B': self.B, 'E': self.E, 'M': self.M, 'correlation_time': self.correlation_time}
        self.observables.append(row)
