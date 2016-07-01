# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:48:40 2016

@author: emc
"""

from solar_extraction import *

dt = 1
#============================================================================#
# Generate SW2 data
# Background = 0
# Constant amplitude and phase
#============================================================================#
data1S = generate_tides('2016-01-01', '2016-01-10', amps=[0,10,10], phase='C',
                        dt=dt, nRange=[2], sRange=[2], filename='SW2_scenario1.txt', 
                        component='solar')

#============================================================================#
# Generate M2 data
# Background = 0
# Constant amplitude and phase
#============================================================================#
data1L = generate_tides('2016-01-01', '2016-01-10', amps=[0,10,10], phase='C',
                        dt=dt, nRange=[2], sRange=[2], filename='M2_scenario1.txt', 
                        component='lunar')

#============================================================================#
# Generate total data
# Background = 0
# Constant amplitude and phase
#============================================================================#
data1T = generate_tides('2016-01-01', '2016-01-10', amps=[0,10,10], phase='C',
                        dt=dt, nRange=[2], sRange=[2], filename='TT_scenario1.txt', 
                        component='s+l')

# Bin by solar local time
means1T = bin_by_solar(data1T)

# Subtract the averages according to solar local time
subtracted, nosol1T = remove_solar(data1T, means1T)

plot_vs_date(nosol1T, -120, 'Lunar tides, dt = {} minutes,'.format(int(dt*60)), 
             data2=data1L, c=['blue', 'deepskyblue'], lb=['reconstructed', 'original'],
             mode='both')


# Bin the results by lunar local time
#result = bin_by_lunar(nosol1T, '1T')
#print('Done with binning by lunar')
# The following plot is done WITHOUT binning by LLT just to see how it looks.
#compare_with_plot(result[:,6], data1L[:,6])
#print('Done')