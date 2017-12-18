import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


def main():
    df = pd.read_csv('test.csv')
    temperature = df['T'].unique()
    N = 25
    energy = np.zeros(len(temperature))
    mag = np.zeros(len(temperature))
    cv = np.zeros(len(temperature))
    chi = np.zeros(len(temperature))

    for i in range(len(temperature)):
        energies = df[df['T'] == temperature[i]]['E'].values
        mags = np.absolute(df[df['T'] == temperature[i]]['M'].values)
        energy[i] = np.average(energies)
        mag[i] = np.average(mags)
        cv[i] = np.var(energies)/(N*temperature[i]**2)
        chi[i] = np.var(mags)/(N*temperature[i])

    pl.plot(temperature, energy)
    pl.show()

    pl.plot(temperature, mag)
    pl.show()
    
    pl.plot(temperature, cv)
    pl.show()

    pl.plot(temperature, chi)
    pl.show()



    
    


if __name__=='__main__':
    main()
