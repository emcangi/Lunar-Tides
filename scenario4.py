# -*- coding: utf-8 -*-
"""
Generate tidal data, bin and average by solar local time, subtract SLT average
from original and bin residual by lunar local time.

Scenario 4:
    Constant background tidal amplitude (user designated)
    Constant SW2 amplitude and phase (user designated)
    Varying M2 amplitude (user designated)
    Varying M2 phase (user designated)

Output (In both cases, {} gets filled in with the time step):
    genM2_dt={}_b{}_sc4.txt         Original generated tidal data
    reconM2_dt={}_b{}_sc4.txt       Tidal data binned by lunar local time and
                                longitude

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from solar_extraction import *
import warnings

# LISTS TO STORE ARRAYS FOR PLOTTING ===========================================
gen_lunar = []                # Arrays of generated lunar tidal data (M2)
recon_lunar = []              # Reconstructed lunar tides after calculations

# PARAMETERS & VARIABLES =======================================================
# You may change these variables
dts = [0.5, 1]               # Time steps for data generation in hours
bin_sz = 0.5                    # Bin size to use when doing SLT and LLT binning
a = [2, 10, 10]               # Tidal amplitudes, format [background, sun, moon]
af = 'M'                      # Amplitude flag
start = '2016-01-01'          # Start date for data generation
end = '2016-01-30'            # End date for data generation
f1 = 'genM2_dt={}_b{}_sc4.txt'    # File to write generated tidal data
f2 = 'reconM2_dt={}_b{}_sc4.txt'  # File to write reconstructed lunar tidal data
L = -180                      # Longitude to use for plotting results


# Do not change the following variables
PHI = 'VM'                    # Constant phase
N = [2]                       # Values of n to use in format [n1, n2...]
S = [2]                       # Values of s to use in format [s1, s2...]


# DO ANALYSIS FOR ALL DT =======================================================
for dt in dts:
    if bin_sz < dt:
        warnings.warn('Warning! Bin size {} < time step {}. This will not '
                      'cause problems with plotting but will cause problems '
                      'with analysis and fitting because some entries in '
                      'either the SLT average or LLT average will be '
                      'zero.'.format(bin_sz, dt))

    # Generate lunar data only (for comparison)
    dataM = generate_tides(start, end, amps=a, ampflag=af, phase=PHI, dt=dt,
                           nrange=N, srange=S, filename=f1.format(dt, bin_sz),
                           component='lunar')

    # Generate total data (for calculation)
    dataT = generate_tides(start, end, amps=a, ampflag=af, phase=PHI, dt=dt,
                           nrange=N, srange=S, component='s+l')

    # append a line including bin size for clarity
    with open(f1.format(dt, bin_sz), 'a') as file:
        file.write('\n bin size = {}\n'.format(bin_sz))

    # Bin generated data by solar local time
    means_slt = bin_by_solar(dataT, bin_sz)

    # Subtract the means by SLT from original data
    nosol = remove_solar(dataT, means_slt, bin_sz)

    # Bin the results by lunar local time
    means_llt = bin_by_lunar(nosol, bin_sz)

    # Write the results of lunar binning to a file
    cells = '{:<20}\t'*3
    line0 = cells.format('Lunar local time', 'Longitude',
                         'Lunar Tide (Reconstructed)')
    np.savetxt(f2.format(dt, bin_sz), means_llt, fmt='%-20.4f', delimiter='\t',
               header=line0, comments='')

    # append a line including bin size for clarity
    with open(f2.format(dt, bin_sz), 'a') as file:
        file.write('\n bin size = {}\n'.format(bin_sz))

    # Insert lunar averages into original data structure (with dates, etc)
    reconM = insert_llt_avgs(nosol, means_llt, bin_sz)

    # build lists of result arrays to use for plotting
    gen_lunar.append(dataM)
    recon_lunar.append(reconM)

# PLOT =========================================================================
# Compares the generated M2 data with the reconstructed M2 data (which uses
# average by LLT)
t = 'Original and reconstructed M2, full lunar cycle, bin={} min, sc4'.format(
    int(60*bin_sz))
plot_vs_date_multi(recon_lunar, L, title=t, dts=dts, data2=gen_lunar,
                   c=['blue', 'deepskyblue'], lb=['Reconstructed M2',
                                                  'Original M2'], mode='both')