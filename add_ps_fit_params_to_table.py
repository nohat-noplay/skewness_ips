import os, math
import numpy as np

from scipy.signal import welch
from scipy.stats import linregress

from astropy.table import Table

from imstack import ImageStack

from iminuit import Minuit
from iminuit.cost import LeastSquares

import argparse

from tqdm import tqdm

import warnings
warnings.filterwarnings('ignore')

def find_ts(table_row, imstack):
    # Returns the time series for a single observation. Table should be the single row for that observation 

    ra = table_row['ra_var'] 
    dec = table_row['dec_var']

    x, y = imstack.world2pix(ra, dec)
    time_series = imstack.pix2ts(x, y, correct=False)

    return time_series 

def get_ps_frac_error(ndof):
    #Returns the fractional error on each power spectrum point
    return np.sqrt(2*ndof)/(ndof*(1 - 2/(9*ndof))**3)

def find_ps(table, imstack):
    #returns the power spectrum, frequency, and errors for one observation

    time_series = find_ts(table,imstack)

    #Remove any overall trend from the time series:
    result = linregress(range(len(time_series)),time_series)
    time_series = np.array(time_series)
    time_series -= result.slope * np.array(range(len(time_series))) + result.intercept

    freq, power_spectrum = welch(time_series, fs=2, window='hanning', nperseg=32, detrend='constant')

    ndof = len(time_series) / len(power_spectrum)  #number of degrees of freedom
    error = get_ps_frac_error(ndof) * (power_spectrum)

    #remove noise (using relationship between off-source noise and moment2_background):
    power_spectrum = power_spectrum - 2.59 * table['moment2_background']**2

    #Remove the first and last values from each
    freq = freq[1:-1] 
    power_spectrum = power_spectrum[1:-1]
    error = error[1:-1]

    return freq, power_spectrum, error

def butterworth_ps(f, knee, alpha, amplitude):
  #model function to fit the power spectrum to
  return amplitude / (1+(f/knee)**(2*alpha))

def find_ps_fit(table, imstack):
    # uses iminuit to fit the model to the power spectrum

    freq, power_spectrum, error = find_ps(table, imstack)

    least_squares = LeastSquares(freq, power_spectrum, error, butterworth_ps)

    m = Minuit(least_squares, knee=table['first_moment'], alpha=2, amplitude=np.mean(power_spectrum[0:len(power_spectrum)//3]))
    m.migrad()
    #m.hesse()

    return m


parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Input file")  
parser.add_argument("outfile", help="Output file") 
args = parser.parse_args()

table = Table.read(args.infile)

#Initialize new columns
for name in ['knee','alpha','amplitude']:
    table[name] = np.nan*np.ones(len(table), dtype=np.float32)
    table[f'{name}_error_lower'] = np.nan*np.ones(len(table), dtype=np.float32)
    table[f'{name}_error_upper'] = np.nan*np.ones(len(table), dtype=np.float32)

table['good_fit'] = np.zeros(len(table), dtype=bool)
table['chi_squared'] = np.nan*np.ones(len(table), dtype=np.float32)

#Count number of bad fits 
failed = 0
large_error = 0

for i, table_row in tqdm(enumerate(table), total=len(table)):
    source_path = os.path.join(table_row['dir'], 
                               str(table_row['obsid'])+".hdf5")
    source_imstack = ImageStack(source_path, freq="121-132")

    try:
        m = find_ps_fit(table_row, source_imstack)
    except Exception as error:
        print(f'error:{error}')
        table['good_fit'][i] = False
        continue

    # check if the fit worked
    if not m.fmin.is_valid:
        print('fit not valid')
        failed += 1
        table['good_fit'][i] = False  
        continue
    
    #Get the errors on the parameters:
    m.minos()
    
    #See whether the errors on the parameters are too large (>50%):
    good = True
    for p in m.params:
        if math.isnan(p.merror[0]) or math.isnan(p.merror[1]) or abs(p.merror[0] / p.value) > 0.5 or abs(p.merror[1] / p.value) > 0.5:
            good = False
            large_error += 1
            break
    
    #enter values into table (even if the errors are >50%)
    for name in m.parameters:
        table[name][i] = m.values[name]
        table[f'{name}_error_lower'][i] = m.merrors[0].lower
        table[f'{name}_error_upper'][i] = m.merrors[0].upper
    table['chi_squared'][i] = m.fmin.reduced_chi2

    #record whether the errors are too large
    table['good_fit'][i] = good

print(f"{failed} fits failed")
print(f"{large_error} errors > 50%")


table.write(args.outfile, format='votable', overwrite=True)
   
