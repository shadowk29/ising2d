from ising2d import ising2d
import numpy as np
import pandas as pd
import progressbar
import itertools


def main():
    temperatures = np.linspace(2.3, 4, 20)
    fields= [0]
    sizes = [40]
    microstates = 1000
    magnet = ising2d(temperatures, fields, sizes, microstates, algorithm = 'wolff', output_folder = 'output', verbose=True)
    magnet.run()

if __name__=='__main__':
    main()
