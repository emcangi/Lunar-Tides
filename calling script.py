# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:48:40 2016

@author: emc
"""

from solar_extraction import *

subs = []
totals = []
dts = [0.1, 0.2, 0.3]

for dt in dts:

    #============================================================================#
    # Generate SW2 data
    # Background = 0sd
    # Constant amplitude and phase
    #============================================================================#
    data1S, SW2only120 = generate_tides('2016-01-01', '2016-01-03',
                                        amps=[0,10,10], phase='C', dt=dt,
                                        nRange=[2], sRange=[2],
                                        filename='SW2_revised.txt',
                                        component='solar')

    #============================================================================#
    # Generate M2 data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    data1M, M2only120 = generate_tides('2016-01-01', '2016-01-03',
                                       amps=[0,10,10], phase='C', dt=dt,
                                       nRange=[2], sRange=[2],
                                       filename='M2_revised.txt',
                                       component='lunar')

    #============================================================================#
    # Generate total data
    # Background = 0
    # Constant amplitude and phase
    #============================================================================#
    data1T, TTonly120 = generate_tides('2016-01-01', '2016-01-03',
                                       amps=[0,10,10], phase='C', dt=dt,
                                       nRange=[2], sRange=[2],
                                       filename='TT_scenario1.txt',
                                        component='s+l')

    # Save a file for one longitude only
    np.savetxt('SW2_revised_lon=-120.txt', SW2only120, fmt='%20.4f',
               delimiter='\t')
    np.savetxt('M2_revised_lon=-120.txt', M2only120, fmt='%20.4f',
               delimiter='\t')
    np.savetxt('TT_revised_lon=-120.txt', TTonly120, fmt='%20.4f',
               delimiter='\t')

    means1T = bin_by_solar(data1T)
    # np.savetxt('means1L_dt={}_revised.txt'.format(dt), means1L, fmt='%-20.4f',
    #            delimiter='\t')


    # Subtract the averages according to solar local time
    avgs, nosol1T = remove_solar(data1T, means1T)

    subs.append(avgs)
    totals.append(data1T)

plot_vs_date_multi(subs, -120, title='Tides and average over SLT, small dt,',
                   dts=dts, data2=totals, c=['red', 'darkorchid'],
                   lb=['Avg by SLT', 'SW2 + M2 tides'],
                   mode='both')

# Bin the results by lunar local time
#result = bin_by_lunar(nosol1T, '1T')
#print('Done with binning by lunar')
# The following plot is done WITHOUT binning by LLT just to see how it looks.
#compare_with_plot(result[:,6], data1L[:,6])
#print('Done')