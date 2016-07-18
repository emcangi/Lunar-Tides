# -*- coding: utf-8 -*-
"""
ANALYZE TIME-GCM DATA

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from lunar_tide_extraction import *
from scipy.stats import chisquare

fname = 'Analysis/TIME-GCM/TUVZ_allLev_timegcm_wg5_zmPlev05_GPI_lunarM2_tmp4.nc'

timegcm_data = wrangle_timegcm_data(fname, 20, 'TN')


### FIX EVERYTHING BELOW HERE !!!!! ----------------------#####################




# =========================== PARAMETERS & VARIABLES ===========================
# You may change these variables
a = [0, 10, 10]                       # Amplitudes: [background, sun, moon]
start = '2016-01-01'                  # Start date for data generation
ends = ['2016-01-15',                 # End date for data generation
        '2016-01-23',
        '2016-01-30']
f1 = 'genM2_{}_dt={}_b{}_sc1.txt'     # File to write generated tidal data
f2 = 'reconM2_{}_dt={}_b{}_sc1.txt'   # File to write reconstructed tidal data
f3 = 'origM2_{}_dt={}_b{}_sc1.txt'    # File to write original M2 data after bin
f4 = 'scenario1_results.txt'          # Prints results of χ² test, etc
L = -120                              # Longitude to use for plotting results
cycle = ['Half', '75%', 'Full']

# Strongly suggested not to change these variables
DTS = [0.5, 1]                # Time steps for data generation in hours
BINS = [0.5, 1]               # Bin size to use when doing SLT and LLT binning
N = [2]                       # Values of v to use in format [n1, n2...]
S = [2]                       # Values of s to use in format [s1, s2...]
GEN_LUNAR = []                # Arrays of generated lunar tidal data (M2)
RECON_LUNAR = []              # Reconstructed lunar tides after calculations


# =================================== MAIN =====================================

# Initialize results file
f = open(f4, 'w')
f.write('SCENARIO 1\n')
f.write('Background tides: {}\n'.format(a[0]))
f.write('Solar amplitude: {}\n'.format(a[1]))
f.write('Solar phase: 0\n')
f.write('Lunar amplitude: {}\n'.format(a[2]))
f.write('Lunar phase: 0\n')
f.write('\n')
f.close()

for end, cyc in zip(ends, cycle):
    for dt in DTS:
        for binsz in BINS:
            # GENERATE DATA ----------------------------------------------------
            if binsz < dt:
                continue


            # Bin generated data by solar local time
            means_slt = bin_by_solar(dataT, binsz)

            # Subtract the means by SLT from original data
            nosol = remove_solar(dataT, means_slt, binsz)

            # Bin the results by lunar local time
            recon_M2_llt_bin = bin_by_lunar(nosol, binsz)

            # Bin original data by LLT -- to compare to reconstruction
            orig_M2_llt_bin = bin_by_lunar(dataM, binsz)

            # Write the results of lunar binning reconstructed M2 to a file
            cells = '{:<20}\t'*3
            line0 = cells.format('Lunar local time', 'Longitude',
                                  'Lunar Tide (Reconstructed)')
            # np.savetxt(f2.format(cyc, dt, binsz), recon_M2_llt_bin,
            # fmt='%-20.4f', delimiter='\t', header=line0, comments='')

            # Write the results of lunar binning original M2 to a file
            # np.savetxt(f3.format(cyc, dt, binsz), orig_M2_llt_bin,
            #            fmt='%-20.4f',
            #            delimiter='\t', header=line0, comments='')

            # χ² ANALYSIS ------------------------------------------------------

            # uses all data across all longitudes
            obs = recon_M2_llt_bin[:, 2]
            exp = orig_M2_llt_bin[:, 2]
            v = len(obs)
            chi, p = chisquare(obs, exp)

            # SCIPY CURVE_FIT --------------------------------------------------
            guess = [a[2], 0]                  # Initial parameter guess
            bounds = [[9, 0], [11, 3*pi/4]]  # [[lo A, lo φ], [hi A, hi φ]]
            ap = amp_and_phase(recon_M2_llt_bin, guess, bounds, N[0], S[0])
            error_amp = round((abs(ap[0] - a[2]) / a[2]) * 100, 2)
            diff_phase = round(ap[1] - 0, 6)

            # WRITE RESULTS ----------------------------------------------------

            with open(f4, 'a') as f:
                f.write('{} lunar cycle\n'.format(cyc))
                f.write('Timestep: {} hr\n'.format(dt))
                f.write('Bin size: {} hr\n'.format(binsz))
                f.write('χ² = {}, p = {} (all longitudes)\n'.format(chi, p))
                f.write('number of data points, v = {}\n'.format(v))
                f.write('χ² ≤ v: {}\n'.format(chi <= v))
                f.write('Amplitudes and phases\n')
                f.write('Average reconstructed lunar amplitude across '
                        'longitudes: {}\n'.format(round(ap[0], 2)))
                f.write('M2 amplitude percent error: {}%\n'.format(error_amp))
                f.write('Average reconstructed lunar phase across longitudes: {'
                        '}\n'.format(round(ap[1], 2)))
                f.write('M2 phase difference from original (0): {}\n'.format(
                    diff_phase))
                f.write('\n')
            f.close()

            # MAKE COMPARISON PLOT ---------------------------------------------

            orig_L = orig_M2_llt_bin[np.where(orig_M2_llt_bin[:, 1] == L)]
            recon_L = recon_M2_llt_bin[np.where(recon_M2_llt_bin[:, 1] == L)]

            # Generate data to plot the fit line
            fit = ap[0] * np.cos((2*pi*N[0] / 24) * recon_L[:, 0] +
                                       (S[0] - N[0]) * L - ap[1])

            plt.figure(figsize=(10,8))
            plt.plot(orig_L[:, 0], orig_L[:, 2], color='deepskyblue',
                     marker='o', markersize=8, label='Original')
            plt.plot(recon_L[:, 0], recon_L[:, 2], color='blue',
                     marker='x', markersize=10, label='Reconstructed')
            plt.plot(recon_L[:, 0], fit, color='red', label='Fit line')
            title = 'M2 vs LLT, {}° longitude, {} cycle, dt={} hr, ' \
                    'b={} hr'.format(L, cyc, dt, binsz)
            plt.title(title)
            plt.xlabel('Lunar local time (hours)')
            plt.ylabel('Tidal amplitude')
            plt.legend(loc='lower right')
            plt.rcParams.update({'font.size': 16})
            plt.tight_layout()
            fn = title + '.png'
            plt.savefig(fn, bbox_inches='tight')
