from ising2d import ising2d
import numpy as np
import pandas as pd


def main():
    temperature = np.linspace(1, 4, 20)
    B = 0
    L = 10
    magnet = ising2d('metropolis')
    for T in temperature:
        magnet.update_system(L, T, B)
        magnet.thermalize()
        magnet.correlation_time()
        for i in range(10000):
            print i
            magnet.update_microstate()
    magnet.print_observables('test.csv')



    
    


if __name__=='__main__':
    main()
