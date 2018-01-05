from ising2d import ising2d
import numpy as np
import pandas as pd
import itertools


def main():
    temperatures = np.linspace(2,4,20)
    fields= [0]
    sizes = [16]
    microstates = 1000
    magnet = ising2d(temperatures, fields, sizes, microstates, algorithm = 'metropolis', output_folder = 'output_test', save_states=10)
    magnet.run()

if __name__=='__main__':
    main()
