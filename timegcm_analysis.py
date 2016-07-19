# -*- coding: utf-8 -*-
"""
ANALYZE TIME-GCM DATA

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from lunar_tide_extraction import *

fname = 'Analysis/TIME-GCM/TUVZ_allLev_timegcm_wg5_zmPlev05_GPI_lunarM2_tmp4.nc'
lat = -45
type = 'temp'

# Extract and format the data from the NetCDF file
timegcm_data = format_timegcm_data(fname, lat, 'TN')
print('Finished formatting data')

# =========================== PARAMETERS & VARIABLES ===========================
# You may change these variables
start = timegcm_data[0, 3]
duration = [15, 23, 30]
fr = 'reconM2_{}_lat={}_timegcm_{}.txt'       # File to write reconstructed data
result_file = 'timegcm_results_{}.txt'        # Stores results
L = 150                                  # Longitude to use for plotting results
cycle = ['Half', '75%', 'Full']

# Strongly suggested not to change these variables
binsz = 1                     # Bin size to use when doing SLT and LLT binning
N = [2]                       # Values of v to use in format [n1, n2...]
S = [2]                       # Values of s to use in format [s1, s2...]


# =================================== MAIN =====================================

# Initialize results file
f = open(result_file.format(type), 'w')
f.write('TIME-GCM Analysis, {}, lat = {}\n'.format(type, lat))
f.close()

for dur, cyc in zip(duration, cycle):

    # Extract only entries that are within the cycle length specified by dur
    data = timegcm_data[np.where(timegcm_data[:, 3] <= start + dur)]

    # Bin generated data by solar local time
    means_slt = bin_by_solar(data, binsz)

    # Subtract the means by SLT from original data
    nosol = remove_solar(data, means_slt, binsz)

    # Bin the results by lunar local time
    recon_M2_llt_bin = bin_by_lunar(nosol, binsz)

    # Write the results of lunar binning reconstructed M2 to a file
    cells = '{:<20}\t'*3
    line0 = cells.format('Lunar local time', 'Longitude',
                         'Lunar Tide (m/s)')
    np.savetxt(fr.format(cyc, lat, type), recon_M2_llt_bin, fmt='%-20.4f',
               delimiter='\t', header=line0, comments='')

    # SCIPY CURVE_FIT --------------------------------------------------
    guess_mean = np.mean(recon_M2_llt_bin[:, 2])
    guess_amp = np.max(recon_M2_llt_bin[:, 2]) - guess_mean
    guess_phase = 0

    guess = [guess_amp, guess_phase, guess_mean]       # Initial parameter guess
    ap = amp_and_phase(recon_M2_llt_bin, guess, N[0], S[0])

    # WRITE RESULTS ----------------------------------------------------

    with open(result_file.format(type), 'a') as f:
        f.write('{} lunar cycle\n'.format(cyc))
        f.write('Bin size: {} hr\n'.format(binsz))
        f.write('Reconstructed lunar amplitude: {}\n'.format(round(ap[0], 2)))
        f.write('Reconstructed lunar phase: {}\n'.format(round(ap[1], 2)))
        f.write('Reconstructed baseline: {}\n'.format(round(ap[2], 2)))
        f.write('\n')
    f.close()

    # MAKE COMPARISON PLOT ---------------------------------------------

    recon_L = recon_M2_llt_bin[np.where(recon_M2_llt_bin[:, 1] == L)]

    # Generate data to plot the fit line
    fit = ap[0] * np.cos((2 * pi * N[0] / 24) * recon_L[:, 0] +
                             (S[0] - N[0]) * L - ap[1]) + ap[2]

    plt.figure(figsize=(10, 8))

    plt.plot(recon_L[:, 0], recon_L[:, 2], color='blue',
             marker='x', markersize=10, label='Reconstructed M2')
    plt.plot(recon_L[:, 0], fit, color='red', label='Fit line')

    title = 'M2 vs LLT, {}° latitude, {}° longitude, {} cycle'.format(lat, L,
                                                                      cyc)
    plt.title(title)
    plt.xlabel('Lunar local time (hours)')
    plt.ylabel('Tidal amplitude (m/s)')
    plt.legend(loc='lower right')
    plt.rcParams.update({'font.size': 16})
    plt.tight_layout()
    fn = title + '.png'
    plt.savefig(fn, bbox_inches='tight')
    #plt.show()
