# -*- coding: utf-8 -*-
"""
ANALYZE TIME-GCM DATA

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from lunar_tide_extraction import *

fname = 'Analysis/TIME-GCM/TUVZ_allLev_timegcm_wg5_zmPlev05_GPI_lunarM2_tmp4.nc'
type = 'zonal wind'

# =========================== PARAMETERS & VARIABLES ===========================
# You may change these variables
duration = 14.5                                    # Cycle to use = 14.5 days
result_file = 'timegcm_results_{}.txt'             # Stores results

# Strongly suggested not to change these variables
BINSZ = 1                                     # Bin size for SLT and LLT binning
N = [2]                                       # Values of n, format [n1, n2...]
S = [2]                                       # Values of s, format [s1, s2...]


# =================================== MAIN =====================================

results = []          # normally bad form, but okay because data isn't that long

for lat in range(-90, 90, 5):

    # Extract and format the data from the NetCDF file
    timegcm_data, real_lat = format_timegcm_data(fname, lat, 'UN')
    LAST_DATE = timegcm_data[-1, 3]  # Ending date in file

    # Set up sliding window
    window_start = timegcm_data[0, 3]              # Obtain Julian starting date
    window_end = window_start + duration           # Window end date
    window_ctr = (window_start + window_end) / 2   # Date in center of window

    while window_end <= LAST_DATE:
        # Extract only entries that are within the current window
        c = (window_start <= timegcm_data[:, 3]) & (timegcm_data[:,
                                                    3] <= window_end)
        data = timegcm_data[np.where(c)]

        # Bin by SLT, subtract SLT average from original data, then bin by LLT
        means_slt = bin_by_solar(data, BINSZ)
        nosol = remove_solar(data, means_slt, BINSZ)
        recon_M2_llt_bin = bin_by_lunar(nosol, BINSZ)

        # SCIPY CURVE_FIT ------------------------------------------------------
        guess_mean = np.mean(recon_M2_llt_bin[:, 2])
        guess_amp = np.max(recon_M2_llt_bin[:, 2]) - guess_mean
        guess_phase = 0

        guess = [guess_amp, guess_phase, guess_mean]   # Initial parameter guess
        ap = amp_and_phase(recon_M2_llt_bin, guess, N[0], S[0])

        # WRITE RESULTS --------------------------------------------------------

        results.append([ap[0], ap[1], ap[2], window_ctr, real_lat])
        window_start += 1
        window_end += 1
        window_ctr += 1
    print('{} complete'.format(lat))

results = np.asarray(results)
cells = '{:<24}'*4
line0 = cells + '\n\n{:<24}' + cells
line0 = line0.format('TIME-GCM Analysis', type, 'Half lunar cycle, 14.5 days',
                     'Bin size: {} hr'.format(BINSZ), 'Amplitude', 'Phase',
                     'Vert offset', 'Date', 'Lat')

np.savetxt(result_file.format(type), results, fmt='%-20.4f',
                       delimiter='\t', header=line0, comments='')


