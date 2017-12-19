from ising2d import ising2d
import numpy as np
import pandas as pd
import progressbar


def main():
    temperatures = np.linspace(1, 4, 20)
    fields= [0]
    sizes = [10]
    microstates = 10000
    magnet = ising2d(temperatures, fields, sizes, microstates, algorithm = 'wolff', output_folder = 'output_wolff', verbose=True)
    magnet.run()

if __name__=='__main__':
    main()
