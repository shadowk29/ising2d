from ising2d import ising2d
import numpy as np
import pandas as pd
import itertools


def main():
    temperatures = [2.0/np.log(1+np.sqrt(2))]
    fields= np.array([10**n for n in np.linspace(-8, -4, 9)])
    sizes = [128]
    microstates = 10000
    magnet = ising2d(temperatures, fields, sizes, microstates, algorithm = 'wolff', output_folder = 'B_output', verbose=True)
    magnet.run()

if __name__=='__main__':
    main()
