import numpy as np
import matplotlib.pyplot as pl
import pandas as pd

class ising2d():
    def __init__(L, T, B):
        self.T = T
        self.B = B
        self.L = L
        self.N = L**2
        self.correlation_time = 10*N
        self.state = np.random.choice([0,1], size=(L,L))
        self.thermalize()

    def thermalize(self):
        pass

    def update_microstate(self):
        pass

    def correlation_time(self):
        pass

    def print_state(self):
        pass

    def print_observables(self):
        pass

    def __spinflip(self, algorithm = 'metropolis'):
        pass

    def print_autocorrelation(self):
        pass

    
    
