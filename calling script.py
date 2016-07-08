# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:48:40 2016

@author: emc
"""

from solar_extraction import *

subs = []
totals = []
dts = [0.1, 0.5, 1]

for dt in dts:

    #============================================================================#
    # Generate SW2 data
    # Background = 0sd
    # Constant amplitude and phase
    #============================================================================#
    data1S, SW2only120 = generate_tides('2016-01-01', '2016-01-15',
                                        amps=[0,10,10], phase='C', dt=dt,
                                        nRange=[2], sRange=[2],
                                        filename='SW2_revised.txt',
                                        component='solar')

    #============================================================================#
    # Generate M2 data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    data1M, M2only120 = generate_tides('2016-01-01', '2016-01-15',
                                       amps=[0,10,10], phase='C', dt=dt,
                                       nRange=[2], sRange=[2],
                                       filename='M2_revised.txt',
                                       component='lunar')

    #============================================================================#
    # Generate total data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    data1T, TTonly120 = generate_tides('2016-01-01', '2016-01-15',
                                       amps=[0,10,10], phase='C', dt=dt,
                                       nRange=[2], sRange=[2],
                                       filename='TT_revised.txt',
                                        component='s+l')

    # # Save a file for one longitude only
    # np.savetxt('SW2_revised_lon=-120.txt', SW2only120, fmt='%20.4f',
    #            delimiter='\t')
    # np.savetxt('M2_revised_lon=-120.txt', M2only120, fmt='%20.4f',
    #            delimiter='\t')
    # np.savetxt('TT_revised_lon=-120.txt', TTonly120, fmt='%20.4f',
    #            delimiter='\t')

    means1M = bin_by_solar(data1M)
    # np.savetxt('means1T_dt={}_revised.txt'.format(dt), means1T, fmt='%-20.4f',
    #             delimiter='\t')


    # Subtract the averages according to solar local time
    avgs, nosol1M = remove_solar(data1M, means1M)
    # longs = np.around(avgs[:, 2], decimals=4)
    # subs120 = avgs[np.where(longs==-120)]
    # np.savetxt('M_subtractedvals_dt={}_revised.txt'.format(dt), avgs,
    #            fmt='%-20.4f',
    #            delimiter='\t')
    # np.savetxt('M_subtracted_values_long=-120.txt', subs120, fmt='%-20.4f',
    #            delimiter='\t')

    subs.append(avgs)
    totals.append(data1M)

plot_vs_date_multi(subs, -120, title='M2 tide and average over SLT, ',
                   dts=dts, data2=totals, c=['blue', 'deepskyblue'],
                   lb=['M2 tide average by SLT', 'M2 original'],
                   mode='both')

# Bin the results by lunar local time
#result = bin_by_lunar(nosol1T, '1T')
#print('Done with binning by lunar')
# The following plot is done WITHOUT binning by LLT just to see how it looks.
#compare_with_plot(result[:,6], data1L[:,6])
#print('Done')