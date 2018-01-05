from ising2d import ising2d
import numpy as np
import pandas as pd
import itertools


def main():
    tau = np.array([10**n for n in np.linspace(-6, -2, 8)])
    temperatures = 2/np.log(1+np.sqrt(2))*(1+tau)
    fields= [0]
    sizes = [8, 16, 32, 64, 128]
    microstates = 1000
    magnet = ising2d(temperatures, fields, sizes, microstates, algorithm = 'wolff', output_folder = 'newdir', save_states=0)
    magnet.run()

if __name__=='__main__':
    main()
