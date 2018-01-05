from ising2d import ising2d
import numpy as np
import pandas as pd
import itertools


def main():
    temperatures = np.linspace(1,4,20)
    fields= [0]
    sizes = [16]
    microstates = 1000
    magnet = ising2d(temperatures, fields, sizes, microstates, algorithm = 'wolff', output_folder = 'output_test')
    magnet.run()

if __name__=='__main__':
    main()
