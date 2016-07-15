# -*- coding: utf-8 -*-
"""
Generate tidal data, bin and average by solar local time, subtract SLT average
from original and bin residual by lunar local time.

Scenario 2:
    Constant background tidal amplitude
    Constant amplitude (user designated)
    Constant phase (user designated)

Output (In both cases, {} gets filled in with the fraction of cycle, timestep or
 bin size):
    genM2_{}_dt={}_b{}_sc2.txt         Original generated tidal data
    reconM2_{}_dt={}_b{}_sc2.txt       Tidal data binned by lunar local time and
                                       longitude
    origM2_{}_dt={}_b{}_sc2.txt        Original tidal data, binned by lunar
                                       local time, with extraneous columns such
                                       as date, SLT, moon phase removed
    scenario2_results.txt              Summarizes each simulation run with
                                       percentage of lunar cycle, bin size, step
                                       size, χ² test on all data and
                                       reconstructed amplitudes and phase

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from solar_extraction import *
from scipy.stats import chisquare

# =========================== PARAMETERS & VARIABLES ===========================
# You may change these variables
a = [2, 10, 10]                       # Amplitudes: [background, sun, moon]
start = '2016-01-01'                  # Start date for data generation
ends = ['2016-01-15',                 # End date for data generation
        '2016-01-23',
        '2016-01-30']
f1 = 'genM2_{}_dt={}_b{}_sc2.txt'     # File to write generated tidal data
f2 = 'reconM2_{}_dt={}_b{}_sc2.txt'   # File to write reconstructed tidal data
f3 = 'origM2_{}_dt={}_b{}_sc2.txt'    # File to write original M2 data after bin
f4 = 'scenario2_results.txt'          # Prints results of χ² test, etc
L = -120                              # Longitude to use for plotting results
cycle = ['Half', '75%', 'Full']       # Simulation length, units of lunar cycle

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
                continue

            # Generate lunar data only (for comparison)
            dataM = generate_tides(start, end, amps=a, dt=dt, nrange=N,
                                   srange=S, component='lunar')

            # Generate total data (for calculation)
            dataT = generate_tides(start, end, amps=a, dt=dt, nrange=N,
                                   srange=S, component='s+l')

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
            # np.savetxt(f3.format(cyc, dt, binsz), orig_M2_llt_bin, fmt='%-20.4f',
            #            delimiter='\t', header=line0, comments='')

            # χ² ANALYSIS ------------------------------------------------------

            # uses all data across all longitudes
            obs = recon_M2_llt_bin[:, 2]
            exp = orig_M2_llt_bin[:, 2]
            v = len(obs)
            chi, p = chisquare(obs, exp)

            # SCIPY CURVE_FIT --------------------------------------------------
            guess = [a[2], 0]  # Initial parameter guess
            bounds = [[9, 0], [11, 2 * pi]]  # [[low A, low φ], [hi A, hi φ]]
            ap = amp_and_phase(recon_M2_llt_bin, guess, bounds, N[0], S[0])

            # Find some averages
            avg = np.mean(ap, axis=0)
            avg_amp = avg[1]
            avg_phase = avg[2]
            error_amp = round((abs(avg_amp - a[2]) / a[2]) * 100, 2)
            diff_phase = round(avg_phase - 0, 6)

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
                        'longitudes: {}\n'.format(round(avg_amp, 2)))
                f.write('M2 amplitude percent error: {}%\n'.format(error_amp))
                f.write('Average reconstructed lunar phase across longitudes: {'
                        '}\n'.format(round(avg_phase, 2)))
                f.write('M2 phase difference from original (0): {}\n'.format(
                    diff_phase))
                f.write('\n')
            f.close()

            # MAKE COMPARISON PLOT ---------------------------------------------

            orig_L = orig_M2_llt_bin[np.where(orig_M2_llt_bin[:, 1] == L)]
            recon_L = recon_M2_llt_bin[np.where(recon_M2_llt_bin[:, 1] == L)]

            # Generate data to plot the fit line
            this_lon = ap[np.where(ap[:, 0] == L)][0]
            fit = this_lon[1] * np.cos((2*pi*N[0] / 24) * recon_L[:, 0] +
                                       (S[0] - N[0]) * L - this_lon[2])

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


            # EXTRA STUFF FOR PLOTTING BY DATE ---------------------------------

            # # Insert LLT averages into original table (with dates, etc)
            # reconM_full_table = insert_llt_avgs(nosol, recon_M2_llt_bin,
            #                                     binsz)
            #
            # # build lists of result arrays to use for plotting
            # gen_lunar.append(dataM)
            # recon_lunar.append(reconM_full_table)

# ========================= DISPLAY RESULTS IN CONSOLE =========================
# with open(f4, 'r') as f:
#     print(f.read())

# =================================== PLOT =====================================
# Compares the generated M2 data with the reconstructed M2 data (which uses
# average by LLT)

# t = 'Original and reconstructed M2, full lunar cycle, binsz={} min, sc2'.format(
#     int(60*binsz))
# plot_vs_date_multi(recon_lunar, L, title=t, dts=dts, data2=gen_lunar,
#                    c=['blue', 'deepskyblue'], lb=['Reconstructed M2',
#                                                   'Original M2'], mode='both')