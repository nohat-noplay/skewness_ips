README

# Exploring Skewness as a Measure of Interplanatary Scintillation
Code Author: Saf Flatters, Data Science CSIRO Intern from Curtin University

Supervisor: John Morgan, Space and Astronomy team (john.morgan@csiro.au)

Feburary 2025


This code analyses MWA IPS GLEAM Survey radio source timeseries observations by calculating true skew values with a Signal to Noise Ratios (SNR) using linear interpolation for each observation (there are multiple observations per radio source). It produces multiple plots of each radio source to visualize the relationship between measured skew, true skew, SNR, elongation and modulation - with error bars. These plots also have lognormal and rice^2 shape models plotted for comparison against an Ordinary Least Squares Line of Best Fit and what I am calling Barlow's Weighted Non-linear Squares Line of Best fit (which takes into account variability of each observation and asymmetric errors). Each Radio Source has a Barlow's Reduced Chi Square value (using the BWNLS model). 
The code provides two exportable tables (Source_Table and Obs_Table) for further analysis.

--------------------------------------------------------------------

## Imports and Files

Saf's docs:
- Gleam_Skewness.ipynb <- this is the run code for this project
- README 
- environment.yml
- requirements.txt


Data: 
- ips/timeseries_full.hdf5 <- a large set of timeseries from the IPS GLEAM survey (both radio sources and offsources)
        (note: timeseries_full.hdf5 is 105MB and therefore not on the repository... see John Morgan for data)
- ssh_to_local.ipynb <-  how I downloaded the large set of timeseries (on and offsource) from the MWA Image Cube. 


Code from others to extract data:
- ips/bright_source_background.vot
- add_ps_fit_params_to_table.py
- bright_cat_pf.fits
- hdf5_to_tables.ipynb
- imstack.py

- Other:
Powerpoint presentation of project also may be uploaded
