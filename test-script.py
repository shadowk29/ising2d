import ising2d
import numpy as np


def main():
    magnet = ising2d.ising2d(5, 3, 0, 'metropolis')
    magnet.thermalize()


if __name__=='__main__':
    main()
