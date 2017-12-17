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
        self.algorithm = algorithm
        self.correlation_time = None
        self.state = np.random.choice([0,1], size=(L,L))
        self.thermalize()

    def thermalize(self):
        for i in range(N**2):
            self.__spinflip()

    def update_microstate(self):
        pass

    def autocorrelation(self, energy):
        xp = ifftshift((energy - np.average(energy))/np.std(energy))
        n = len(xp)
        xp = np.r_[xp[:n/2], np.zeros_like(xp), xp[n/2:]]
        f = fft(xp)
        S = np.absolute(f)**2
        R = ifft(S)
        self.autocorrelation = np.real(R)[:n/2]/(np.arange(n/2)[::-1]+n/2)
    
    def correlation_time(self, plot=False):
        delays = np.arange(len(self.autocorrelation))
        popt, pcov = curve_fit(exponential, delays, self.autocorrelation, p0=[self.N])
        self.correlation_time = popt[0]
        if plot:

    def print_state(self):
        pass

    def print_observables(self):
        pass

    def __spinflip(self, algorithm = 'metropolis'):
        pass

    def print_autocorrelation(self):
        pass

    def exponential(self, y, n0):
        n = np.arange(len(y))
        return np.exp(-(n/n0))

    
    
