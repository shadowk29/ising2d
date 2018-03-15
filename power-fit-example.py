import numpy as np
import matplotlib.pyplot as pl
from scipy.stats import t as student_t

def linear_power_fit(x, y):
    try:
        ##linearize the model by log-transforming both sides
        popt, pcov = np.polyfit(np.log(x), np.log(y), 1, cov=True)
        power = popt[0]
        prefactor = np.exp(popt[1])
        perr = np.sqrt(np.diag(pcov))[0]
    except ValueError:
        raise ValueError('You need at least 5 data points to get an error estimate for a linear fit')
    return power, prefactor, perr

def generate_data(A, p):
    ##generate some noisy test data
    x = np.arange(1,6)
    y = A*x**p + np.random.normal(0, 1, len(x))
    return x, y

def main():
    x, y = generate_data(5, 2)
    p, A, perr = linear_power_fit(x, y)

    pl.loglog(x, y, 'o', label='Actual power: {0}'.format(2))
    pl.loglog(x, A*x**p, label=r'Fitted power: {0:.2f}$\pm${1:.2f}'.format(p, perr))
    pl.xlabel('X')
    pl.ylabel('Y')
    pl.legend(loc='best')
    pl.tight_layout()
    pl.show()
    
if __name__=='__main__':
    main()
