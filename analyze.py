import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


def main():
    df = pd.read_csv('test.csv')
    temperature = df['T'].unique()
    N = 100
    energy = np.zeros(len(temperature))
    mag = np.zeros(len(temperature))
    cv = np.zeros(len(temperature))
    chi = np.zeros(len(temperature))

    for i in range(len(temperature)):
        energies = df[df['T'] == temperature[i]]['E'].values
        mags = np.absolute(df[df['T'] == temperature[i]]['M'].values)
        energy[i] = np.average(energies)/N
        mag[i] = np.average(mags)/N
        cv[i] = np.var(energies)/(N*temperature[i]**2)
        chi[i] = np.var(mags)/(N*temperature[i])

    pl.plot(temperature, energy, 'o-')
    pl.axvline(x=2/np.log(1+np.sqrt(2)))
    pl.show()

    pl.plot(temperature, mag, 'o-', np.linspace(1,4,200), (1-np.sinh(2.0/np.linspace(1,4,200))**(-4.0))**(1.0/8.0))
    pl.axvline(x=2/np.log(1+np.sqrt(2)))
    pl.show()
    
    pl.plot(temperature, cv, 'o-')
    pl.axvline(x=2/np.log(1+np.sqrt(2)))
    pl.show()

    pl.plot(temperature, chi, 'o-')
    pl.axvline(x=2/np.log(1+np.sqrt(2)))
    pl.show()



    
    


if __name__=='__main__':
    main()
