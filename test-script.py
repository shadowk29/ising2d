import ising2d
import numpy as np


def main():
    magnet = ising2d.ising2d(10, 1.5, 0, 'metropolis')
    magnet.thermalize()
    magnet.correlation_time(plot=True)
    print magnet.corrtime
    for i in range(10):
        magnet.update_microstate()
        print magnet.E


if __name__=='__main__':
    main()
