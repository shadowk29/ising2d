import ising2d
import numpy as np


def main():
    magnet = ising2d.ising2d(4, 2.26918531421302196811, 0, 'metropolis')
    magnet.thermalize()
    magnet.correlation_time(plot=True)
    print magnet.corrtime
    print magnet.probability
    for i in range(100):
        magnet.update_microstate()
        magnet.save_observables()
    magnet.print_observables('test.csv')
    


if __name__=='__main__':
    main()
