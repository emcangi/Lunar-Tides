# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:48:40 2016

@author: emc
"""

from solar_extraction import *

subs = []
totals = []
dts = [0.1]#, 0.25, 1]

for dt in dts:

    #============================================================================#
    # Generate SW2 data
    # Background = 0sd
    # Constant amplitude and phase
    #============================================================================#
    # data1S = generate_tides('2016-01-01', '2016-01-02', amps=[0,10,10], phase='C',
    #                         dt=dt, nRange=[1,2], sRange=[1,2], component='solar')

    #============================================================================#
    # Generate M2 data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    data1L, M2only120 = generate_tides('2016-01-01', '2016-01-05', amps=[0,
                                                                         10,10],
                            phase='C', dt=dt, nRange=[2], sRange=[2],
                            filename='M2_dt={}_UT_eq.txt'.format(dt),
                            component='lunar')

    #============================================================================#
    # Generate total data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    # data1T = generate_tides('2016-01-01', '2016-01-02', amps=[0,10,10], phase='C',
    #                         dt=dt, nRange=[2], sRange=[2], filename='TT_scenario1.txt',
    #                         component='s+l')

    # Bin by local time

    # Save a file
    np.savetxt('M2_lon-120_dt={}_UT_eq.txt'.format(dt), M2only120,
               fmt='%-20.4f', delimiter='\t')

    means1L = bin_by_solar(data1L)
    np.savetxt('means1L_dt={}_UT_eq.txt'.format(dt), means1L, fmt='%-20.4f',
               delimiter='\t')


    # Subtract the averages according to solar local time
    avgs, nosol1T = remove_solar(data1L, means1L)

    subs.append(avgs)
    totals.append(data1L)

plot_vs_date_multi(subs, -120, title='Lunar semidiurnal tide and average over '
                    'SLT,', dts=dts, data2=totals, c=['red', 'darkorchid'],
                   lb=['Lunar avg by SLT', 'lunar semidiurnal original'],
                   mode='show')

# Bin the results by lunar local time
#result = bin_by_lunar(nosol1T, '1T')
#print('Done with binning by lunar')
# The following plot is done WITHOUT binning by LLT just to see how it looks.
#compare_with_plot(result[:,6], data1L[:,6])
#print('Done')