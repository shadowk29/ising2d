from ising2d import ising2d
import numpy as np
import pandas as pd
import progressbar


def main():
    temperature = np.linspace(1, 4, 20)
    B = 0
    L = 10
    states = 10000
    magnet = ising2d('metropolis')
    with progressbar.ProgressBar(max_value=len(temperature)*states) as bar:
        for T in temperature:
            bar.update()
            magnet.update_system(L, T, B)
            magnet.thermalize()
            magnet.correlation_time()
            for i in range(states):
                bar.update(np.where(temperature == T)[0][0]*states + i)
                magnet.update_microstate()
    magnet.print_observables('test.csv')

if __name__=='__main__':
    main()
