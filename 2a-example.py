import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


def main():
    ## read in the data you generated
    observables = pd.read_csv('observables.csv')

    ## get a list of unique temperatures and system sizes, and print them to check
    temperatures = observables['T'].unique()
    lengths = observables['L'].unique()


    ## define critical temperature for use later
    Tc=2/np.log(1+np.sqrt(2))

    ## define an array to hold the energy as a function of temperature, and the standard deviation
    energy = np.zeros(len(temperatures))
    energy_err = np.zeros(len(temperatures))

    ## the main loop. We want to get E(T) for each length
    for L in lengths:
        N = L**2
        ## select the subset of rows of the dataset which have the right L value
        subset = observables[observables['L'] == L]
        ## loop over every temperature
        for i in range(len(temperatures)):
            ## select the energy column for rows of the subset which have the right T value
            ## note that because T is not an integer, we cannot use == to compare.
            T = temperatures[i]
            energies = subset[np.absolute(subset['T'] - T) < 1e-9]['E'].values
            ## be sure to normalize everything by N!
            energy[i] = np.average(energies)/float(N)
            energy_err[i] = np.std(energies)/float(N)
        ##plot the result, one line for each system size
        pl.errorbar(temperatures, energy, yerr=energy_err, fmt='o', label='L={0}'.format(L))

    ## finalize plotting options
    pl.xlabel('Temperature')
    pl.ylabel('Energy per Spin')
    pl.legend(loc='best')
    pl.tight_layout()
    pl.show()
    
    
if __name__=='__main__':
    main()
