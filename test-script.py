from ising2d import ising2d
import numpy as np
import pandas as pd
import itertools


def main():
    temperatures = [2.0/np.log(1+np.sqrt(2))]
    fields= [0.05]
    sizes = [64]
    microstates = 1000
    magnet = ising2d(temperatures, fields, sizes, microstates, algorithm = 'wolff', output_folder = 'output_test')
    magnet.run()

if __name__=='__main__':
    main()
