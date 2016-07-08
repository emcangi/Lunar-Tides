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
    #==========================================================================#
    # Generate SW2 data
    # Background = 0sd
    # Constant amplitude and phase
    #==========================================================================#
    data1S, SW2_120 = generate_tides('2016-01-01', '2016-01-04',
                                     amps=[0, 10, 10], phase='C', dt=dt,
                                     nRange=[2], sRange=[2],
                                     filename='SW2_revised_dt{}.txt'.format(dt),
                                     component='solar')

    #============================================================================#
    # Generate M2 data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    data1M, M2_120 = generate_tides('2016-01-01', '2016-01-04',
                                    amps=[0,10,10], phase='C', dt=dt,
                                    nRange=[2], sRange=[2],
                                    filename='M2_revised_dt{}.txt'.format(dt),
                                    component='lunar')

    #============================================================================#
    # Generate total data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    data1T, TT_120 = generate_tides('2016-01-01', '2016-01-04',
                                    amps=[0,10,10], phase='C', dt=dt,
                                    nRange=[2], sRange=[2],
                                    filename='TT_revised_dt{}.txt'.format(dt),
                                    component='s+l')

    # Save a file for one longitude only
    # np.savetxt('SW2_revised_lon=-120_dt{}.txt'.format(dt), SW2_120,
    #            fmt='%20.4f', delimiter='\t')
    # np.savetxt('M2_revised_lon=-120.txt', M2_120, fmt='%20.4f',
    #            delimiter='\t')
    # np.savetxt('TT_revised_lon=-120.txt', TT_120, fmt='%20.4f',
    #            delimiter='\t')

    means = bin_by_solar(data1S, 0.5)
    # np.savetxt('means_dt={}.txt'.format(dt), means, fmt='%-20.4f',
    #             delimiter='\t')


    # Subtract the averages according to solar local time
    avgs, nosol1 = remove_solar(data1S, means, 0.5)
    # longs = np.around(avgs[:, 2], decimals=4)
    # subs120 = avgs[np.where(longs==-120)]
    # np.savetxt('SW2_lon=-120_dt={}_subtracts_restricted_to_120.txt'.format(dt),
    #            avgs, fmt='%-20.4f', delimiter='\t')
    #
    subs.append(avgs)
    totals.append(data1S)

plot_vs_date_multi(subs, -120, title='SW2 tide and average over SLT, '
                                     'bin=30min,',
                   dts=dts, data2=totals, c=['red', 'orange'],
                   lb=['SW2 tide average by SLT', 'SW2 original'],
                   mode='both')

# Bin the results by lunar local time
#result = bin_by_lunar(nosol1T, '1T')
#print('Done with binning by lunar')
# The following plot is done WITHOUT binning by LLT just to see how it looks.
#compare_with_plot(result[:,6], data1L[:,6])
#print('Done')