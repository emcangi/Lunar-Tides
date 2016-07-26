# -*- coding: utf-8 -*-
"""
ANALYZE TIME-GCM DATA - SLT ONLY

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from lunar_tide_extraction import *

# =========================== PARAMETERS & VARIABLES ===========================
# You may change these variables
duration = 14.5                                    # Cycle to use = 14.5 days
result_file = 'timegcm_results_SLT_{}121.txt'             # Stores results
tide_type = 'zonal wind'
fname = 'Analysis/TIME-GCM/U_plev-4_timegcm_wg5_zmPlev05_lunar.nc'
#fname = 'Analysis/TIME-GCM/TUVZ_allLev_timegcm_wg5_zmPlev05_GPI_lunarM2_tmp4
# .nc'
z = 121

# Strongly suggested not to change these variables
BINSZ = 1                                     # Bin size for SLT and LLT binning
N = [2]                                       # Values of n, format [n1, n2...]
S = [2]                                       # Values of s, format [s1, s2...]

if tide_type == 'temp':
    v = 'TN'
elif tide_type == 'zonal wind':
    v = 'UN'
else:
    raise ValueError("Type must be specified as either 'temp' or 'zonal wind'!")


# =================================== MAIN =====================================

results = []          # normally bad form, but okay because data isn't that long

for lat in range(-90, 90, 5):

    # Extract and format the data from the NetCDF file
    timegcm_data, real_lat = format_timegcm_data(fname, lat, v)
    LAST_DATE = timegcm_data[-1, 3]  # Ending date in file
    START_DATE = timegcm_data[0, 3]

    # Set up sliding window
    window_start = START_DATE              # Obtain Julian starting date
    window_end = window_start + duration           # Window end date
    window_ctr = (window_start + window_end) / 2   # Date in center of window

    while window_end <= LAST_DATE:
        # Extract only entries that are within the current window
        in_window = (window_start <= timegcm_data[:, 3]) & \
                    (timegcm_data[:, 3] <= window_end)
        data = timegcm_data[np.where(in_window)]

        # Bin by SLT, subtract SLT average from original data, then bin by LLT
        means_slt = bin_by_solar(data, BINSZ)

        # SCIPY CURVE_FIT ------------------------------------------------------
        guess_mean = np.mean(means_slt[:, 2])
        max_A = np.max(means_slt[:, 2])
        guess_amp = max_A - guess_mean
        guess_phase = 0
        guess = [guess_amp, guess_phase, guess_mean]   # Initial parameter guess

        bounds = [[0, -pi/2, -np.inf], [15000, pi/2, np.inf]]
        ap = fit_m2(means_slt, guess, N[0], S[0], bounds=bounds)

        # WRITE RESULTS --------------------------------------------------------

        results.append([ap[0], ap[1], ap[2], window_ctr, real_lat])
        window_start += 1
        window_end += 1
        window_ctr += 1

    print('{} complete'.format(lat))

results = np.asarray(results)
cells = '{:<24}'*4
line0 = cells + '\n\n{:<24}' + cells
line0 = line0.format('TIME-GCM Analysis - SOLAR COMPONENT', tide_type,
                     'Half lunar cycle, 14.5 days', 'Bin size: {} '
                     'hr'.format(BINSZ), 'Amplitude', 'Phase', 'Vert offset',
                     'Date', 'Lat')

np.savetxt(result_file.format(tide_type), results, fmt='%-20.4f',
           delimiter='\t', header=line0, comments='')

# # ============= MANIPULATE DATA TO PREPARE FOR PLOTTING ========================
# # Create new array of data in a logical format
# amp_data = np.column_stack((results[:, 3], results[:, 4], results[:, 0]))
# phase_data = np.column_stack((results[:, 3], results[:, 4], results[:, 1]))
#
# # Change numbers so they function as indices
# temp = np.copy(amp_data)
# temp[:, 0] = temp[:, 0] - ((START_DATE + duration) / 2)
# temp[:, 1] = temp[:, 1] + 88.75
#
# temp2 = np.copy(phase_data)
# temp2[:, 0] = temp2[:, 0] - ((START_DATE + duration) / 2)
# temp2[:, 1] = temp2[:, 1] + 88.75
#
# amps = np.zeros([int(max(temp[:, 1]) + 1), int(max(temp[:, 0]) + 1)])
# phases = np.zeros([int(max(temp[:, 1]) + 1), int(max(temp[:, 0]) + 1)])
#
# # Fill the array for plotting with values according to location - this fills
# # it upside down because of the way arrays vs plots are indexed
# for a, p in zip(temp, temp2):
#     x = a[0]
#     y = a[1]
#     z1 = a[2]
#     z2 = p[2]
#     amps[y, x] = z1
#     phases[y, x] = z2
#
# amps = amps[~np.all(amps == 0, axis=1)]  # Get rid of lines with all 0s
# phases = phases[~np.all(phases == 0, axis=1)]
#
# #=============================== MAKE PLOT =====================================
#
# # Make tick labels
# column_labels = np.asarray(sorted(set(amp_data[:, 0])))
# column_labels = np.around(column_labels - START_DATE, 2)
# row_labels = list(np.around(sorted(set(amp_data[:, 1])), 2))
#
# # Blank out most labels so they don't clutter things
# for i in range(len(row_labels)):
#     if i % 5 != 0:
#         row_labels[i] = ''
#
# # Plot a heat map
# fig = plt.figure(figsize=(16, 10))
# plt.axis([0, 41, 0, 36])  # this line makes sure the axes fit the fig!
# heatmap_A = plt.pcolor(amps, cmap=plt.cm.seismic)
#
# ax = fig.gca()
# ax.set_xticks(np.arange(0, amps.shape[1]))
# ax.set_yticks(np.arange(0, amps.shape[0]))
#
# ax.set_xticklabels(column_labels, rotation=45)
# ax.set_yticklabels(row_labels)
# plt.ylabel('Latitude')
# plt.xlabel('Days since 2013-01-01')
#
# plt.rcParams.update({'font.size': 16})
# plt.rc('xtick', labelsize=8)
# plt.rc('ytick', labelsize=12)
# plt.title('SW2 tidal amplitude (zonal wind speeds), z={} km'.format(z), y=1.05)
#
# maxtide = round(np.max(amps))
# cbarticks = np.linspace(0, maxtide, 11)
# cbarticklabels = ['{} m/s'.format(round(i/100, 2)) for i in cbarticks]
# cbar = fig.colorbar(heatmap_A, ticks=cbarticks)
# cbar.ax.set_yticklabels(cbarticklabels)
#
# plt.savefig('sw2tidebylat_amp.png', facecolor='w', edgecolor='w',
#             format='png', transparent=False, bbox_inches='tight')
#
# plt.show()
#
# fig = plt.figure(figsize=(16, 10))
# plt.axis([0, 41, 0, 36])  # this line makes sure the axes fit the fig!
# heatmap_P = plt.pcolor(phases, cmap=plt.cm.seismic)
#
# ax = fig.gca()
# ax.set_xticks(np.arange(0, phases.shape[1]))
# ax.set_yticks(np.arange(0, phases.shape[0]))
#
# ax.set_xticklabels(column_labels, rotation=45)
# ax.set_yticklabels(row_labels)
# plt.ylabel('Latitude')
# plt.xlabel('Days since 2013-01-01')
#
# plt.rcParams.update({'font.size': 16})
# plt.rc('xtick', labelsize=8)
# plt.rc('ytick', labelsize=12)
# plt.title('SW2 tidal phases (zonal wind speeds), z={} km'.format(z), y=1.05)
#
# maxtide = round(np.max(phases))
# mintide = np.min(phases)
# cbarticks = np.linspace(mintide, maxtide, 10)
# cbarticklabels = ['{}°'.format(round(i*180/pi, 2)) for i in cbarticks]
# cbar = fig.colorbar(heatmap_P, ticks=cbarticks)
# cbar.ax.set_yticklabels(cbarticklabels)
#
# plt.savefig('sw2tidebylat_phase.png', facecolor='w', edgecolor='w',
#             format='png', transparent=False, bbox_inches='tight')
#
# plt.show()
