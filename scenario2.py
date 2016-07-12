# -*- coding: utf-8 -*-
"""
Generate tidal data, binsz and average by solar local time, subtract SLT average
from original and binsz residual by lunar local time.

Scenario 2:
    Varying background tidal amplitude
    Constant amplitude (user designated)
    Constant phase (user designated)

Output (In both cases, {} gets filled in with the time step or binsz size):
    genM2_dt={}_b{}_sc2.txt         Original generated tidal data
    reconM2_dt={}_b{}_sc2.txt       Tidal data binned by lunar local time and
                                    longitude

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from solar_extraction import *
import warnings
from scipy.stats import chisquare

# =========================== PARAMETERS & VARIABLES ===========================
# You may change these variables
a = [2, 10, 10]                       # Amplitudes: [background, sun, moon]
af = 'B'                              # Amplitude flag
start = '2016-01-01'                  # Start date for data generation
ends = ['2016-01-15',                 # End date for data generation
        '2016-01-23',
        '2016-01-30']
f1 = 'genM2_{}_dt={}_b{}_sc2.txt'     # File to write generated tidal data
f2 = 'reconM2_{}_dt={}_b{}_sc2.txt'   # File to write reconstructed tidal data
f3 = 'origM2_{}_dt={}_b{}_sc2.txt'    # File to write original M2 data after binsz
f4 = 'scenario2_results.txt'          # Prints results of χ² test, etc
L = -120                              # Longitude to use for plotting results
cycle = ['Half', '75%', 'Full']       # Simulation length, units of lunar cycle

# Strongly suggested not to change these variables
DTS = [0.5, 1]                # Time steps for data generation in hours
BINS = [0.5, 1]             # Bin size to use when doing SLT and LLT binning
PHI = 'C'                     # Constant phase
N = [2]                       # Values of n to use in format [n1, n2...]
S = [2]                       # Values of s to use in format [s1, s2...]
GEN_LUNAR = []                # Arrays of generated lunar tidal data (M2)
RECON_LUNAR = []              # Reconstructed lunar tides after calculations


# =================================== MAIN =====================================

# Initialize results file
f = open(f4, 'w')
f.write('SCENARIO 2\n')
f.write('Background tides: {} (maximum; actual tide varies '
        'sinusoidally)\n'.format(a[0]))
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
            dataM = generate_tides(start, end, amps=a, ampflag=af, phase=PHI,
                                   dt=dt, nrange=N, srange=S, component='lunar')

            # Generate total data (for calculation)
            dataT = generate_tides(start, end, amps=a, ampflag=af, phase=PHI,
                                   dt=dt, nrange=N, srange=S, component='s+l')

            # Bin generated data by solar local time
            means_slt = bin_by_solar(dataT, binsz)

            # Subtract the means by SLT from original data
            nosol = remove_solar(dataT, means_slt, binsz)

            # Bin the results by lunar local time
            recon_M2_llt_bin = bin_by_lunar(nosol, binsz)

            # Bin original data by LLT -- to compare to reconstruction
            orig_M2_llt_bin = bin_by_lunar(dataM, binsz)

            # Write the results of lunar binning reconstructed M2 to a file
            cells = '{:<20}\t' * 3
            line0 = cells.format('Lunar local time', 'Longitude',
                                 'Lunar Tide (Reconstructed)')
            # np.savetxt(f2.format(cyc, dt, bin_sz), recon_M2_llt_bin,
            # fmt='%-20.4f',
            #             delimiter='\t', header=line0, comments='')

            # Write the results of lunar binning original M2 to a file
            np.savetxt(f3.format(cyc, dt, binsz), orig_M2_llt_bin, fmt='%-20.4f',
                       delimiter='\t', header=line0, comments='')

            # χ² ANALYSIS ------------------------------------------------------

            # calculate χ²
            obs = recon_M2_llt_bin[:, 2]
            exp = orig_M2_llt_bin[:, 2]
            chi, p = chisquare(obs, exp)

            # WRITE RESULTS ----------------------------------------------------

            with open(f4, 'a') as f:
                if flag:
                    f.write('WARNING! binsz size is smaller than timestep. '
                            'Analysis may be compromised.\n')
                f.write('{} lunar cycle\n'.format(cyc))
                f.write('Timestep: {} hr\n'.format(dt))
                f.write('Bin size: {}\n'.format(binsz))
                f.write('χ² test: χ² = {}, p = {}\n'.format(chi, p))
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

# t = 'Original and reconstructed M2, full lunar cycle, binsz={} min, sc2'.format(
#     int(60*binsz))
# plot_vs_date_multi(recon_lunar, L, title=t, dts=dts, data2=gen_lunar,
#                    c=['blue', 'deepskyblue'], lb=['Reconstructed M2',
#                                                   'Original M2'], mode='both')