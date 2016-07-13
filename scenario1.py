# -*- coding: utf-8 -*-
"""
Generate tidal data, binsz and average by solar local time, subtract SLT average
from original and binsz residual by lunar local time.

Scenario 1:
    No background tidal amplitudes
    Constant amplitude (user designated)
    Constant phase (zero)

Output (In both cases, {} gets filled in with the time step):
    genM2_dt={}_sc1.txt         Original generated tidal data
    reconM2_dt={}_sc1.txt       Tidal data binned by lunar local time and
                                longitude

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from solar_extraction import *
import warnings
from scipy.stats import chisquare
from scipy.optimize import curve_fit


# =========================== PARAMETERS & VARIABLES ===========================
# You may change these variables
a = [0, 10, 10]                       # Amplitudes: [background, sun, moon]
start = '2016-01-01'                  # Start date for data generation
ends = ['2016-01-15',                 # End date for data generation
       '2016-01-23',
       '2016-01-30']
f1 = 'genM2_{}_dt={}_b{}_sc1.txt'     # File to write generated tidal data
f2 = 'reconM2_{}_dt={}_b{}_sc1.txt'   # File to write reconstructed tidal data
f3 = 'origM2_{}_dt={}_b{}_sc1.txt'    # File to write original M2 data after binsz
f4 = 'scenario1_results.txt'          # Prints results of χ² test, etc
L = -120                              # Longitude to use for plotting results
cycle = ['Half', '75%', 'Full']

# Strongly suggested not to change these variables
DTS = [0.5, 1]                # Time steps for data generation in hours
BINS = [0.5, 1]               # Bin size to use when doing SLT and LLT binning
PHI = 'C'                     # Constant phase
N = [2]                       # Values of n to use in format [n1, n2...]
S = [2]                       # Values of s to use in format [s1, s2...]
GEN_LUNAR = []                # Arrays of generated lunar tidal data (M2)
RECON_LUNAR = []              # Reconstructed lunar tides after calculations


# =================================== MAIN =====================================

# Initialize results file
f = open(f4, 'w')
f.write('SCENARIO 1\n')
f.write('Background tides: 0\n')
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
                warnings.warn('Warning! Bin size {} < time step {}. This will '
                              'not cause problems with plotting but will cause '
                              'problems with analysis and fitting because some '
                              'entries in either the SLT average or LLT '
                              'average will be zero.'.format(binsz, dt))
                flag = True
            else:
                flag = False

            # Generate lunar data only (for comparison)
            dataM = generate_tides(start, end, amps=a, phase=PHI, dt=dt,
                                   nrange=N, srange=S, component='lunar')

            # Generate total data (for calculation)
            dataT = generate_tides(start, end, amps=a, phase=PHI, dt=dt,
                                   nrange=N, srange=S, component='s+l')

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

            # calculate χ² for reconstructed vs. original data across all
            # longitudes
            obs = recon_M2_llt_bin[:, 2]
            exp = orig_M2_llt_bin[:, 2]
            n = len(obs)
            chi, p = chisquare(obs, exp)

            # SCIPY CURVE_FIT --------------------------------------------------
            # Fitting function
            def lunartide(llt, A, P):
                """
                Simple fitting function for the lunar tide.
                ---INPUT---
                :param llt: lunar local time (independent variable)
                :param A: amplitude
                :param P: phase
                :return: array of values
                """
                n = 2
                s = 2
                global L
                return A * np.cos((2 * pi * n / 24) * llt - (n + s) * L + P)

            # Grab data

            time = orig_M2_llt_bin[:]
            recon_M2_llt_bin
            # Initial guess.
            p0 = np.array([10, 0.01])
            L = -180                      # set longitude - should be a loop

            popt, pcov = curve_fit(lunartide, orig180[:, 0], orig180[:, 2], p0,
                                   bounds=([9, 0], [11, 2 * pi]))

            # WRITE RESULTS ----------------------------------------------------

            with open(f4, 'a') as f:
                if flag:
                    f.write('WARNING! binsz size is smaller than timestep. '
                            'Analysis may be compromised.\n')
                f.write('{} lunar cycle\n'.format(cyc))
                f.write('Timestep: {} hr\n'.format(dt))
                f.write('Bin size: {}\n'.format(binsz))
                f.write('χ² test: χ² = {}, p = {}\n'.format(chi, p))
                f.write('n = {}'.format(n))
                f.write('χ² ≤ n: {}'.format(chi <= n))
                f.write('\n')
            f.close()

            # EXTRA STUFF FOR PLOTTING BY DATE ---------------------------------

            # # Insert LLT averages into original table (with dates, etc)
            # reconM_full_table = insert_llt_avgs(nosol, recon_M2_llt_bin,
            #                                     binsz)
            #
            # # build lists of result arrays to use for plotting
            # gen_lunar.append(dataM)
            # recon_lunar.append(reconM_full_table)


# ========================= DISPLAY RESULTS IN CONSOLE =========================
with open(f4, 'r') as f:
    print(f.read())


# =================================== PLOT =====================================
# Compares the generated M2 data with the reconstructed M2 data (which uses
# average by LLT)

# t = 'Original and reconstructed M2, full lunar cycle, binsz={} min, sc1'.format(
#     int(binsz * 60))
# plot_vs_date_multi(recon_lunar, L, title=t, dts=dts, data2=gen_lunar,
#                    c=['blue', 'deepskyblue'], lb=['Reconstructed M2',
#                                                   'Original M2'], mode='show')