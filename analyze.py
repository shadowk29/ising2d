import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


def main():
    df = pd.read_csv('observables.csv')
    temperature = df['T'].unique()
    N = 100
    energy = np.zeros(len(temperature))
    mag = np.zeros(len(temperature))
    cv = np.zeros(len(temperature))
    chi = np.zeros(len(temperature))
    corrtime = np.zeros(len(temperature))
    max_cluster_size = np.zeros(len(temperature))
    cluster_count = np.zeros(len(temperature))
    for i in range(len(temperature)):
        energies = df[df['T'] == temperature[i]]['E'].values
        mags = np.absolute(df[df['T'] == temperature[i]]['M'].values)
        corrtimes = df[df['T'] == temperature[i]]['correlation_time'].values
        clustercounts = df[df['T'] == temperature[i]]['cluster_count'].values
        max_cluster_sizes = df[df['T'] == temperature[i]]['max_cluster_size'].values
        
        energy[i] = np.average(energies)/N
        mag[i] = np.average(mags)/N
        cv[i] = np.var(energies)/(N*temperature[i]**2)
        chi[i] = np.var(mags)/(N*temperature[i])
        corrtime[i] = np.average(corrtimes)
        max_cluster_size[i] = np.average(max_cluster_sizes)
        cluster_count[i] = np.average(clustercounts)

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

    pl.plot(temperature, corrtime, 'o-')
    pl.axvline(x=2/np.log(1+np.sqrt(2)))
    pl.show()

    pl.plot(temperature, cluster_count, 'o-')
    pl.axvline(x=2/np.log(1+np.sqrt(2)))
    pl.show()

    pl.plot(temperature, max_cluster_size, 'o-')
    pl.axvline(x=2/np.log(1+np.sqrt(2)))
    pl.show()


if __name__=='__main__':
    main()
