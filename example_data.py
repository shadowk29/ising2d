from ising2d import ising2d
import numpy as np

def main():
    tauneg = -np.array([10**n for n in np.linspace(-2, -0.26, 8)])
    taupos = np.array([10**n for n in np.linspace(-2, -0.12, 8)])
    tau = np.append(tauneg, taupos)
    tau = np.append(tau, 0)
    temperatures = 2/np.log(1+np.sqrt(2))*(1+tau)
    fields= [0]
    sizes = [32, 8]
    microstates = 10
    magnet = ising2d(temperatures, fields, sizes, microstates, output_folder = 'Example_Data', save_states=10)
    magnet.run()

if __name__=='__main__':
    main()
