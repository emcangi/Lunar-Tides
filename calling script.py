# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:48:40 2016

@author: emc
"""

from solar_extraction import *


#============================================================================#
# Generate SW2 data
# Background = 0
# Constant amplitude and phase
#============================================================================#
data1S = generate_tides('2016-01-01', '2016-03-31', amps=[0,10,10], phase='C',
                        longIncr=15, nRange=[2], sRange=[2], 
                        filename='SW2_0bg_CACP.txt', component='solar')

    
#============================================================================#
# Generate M2 data
# Background = 0
# Constant amplitude and phase
#============================================================================#
data1L = generate_tides('2016-01-01', '2016-03-31', amps=[0,10,10], phase='C',
                        longIncr=15, nRange=[2], sRange=[2], 
                        filename='M2_0bg_CACP.txt', component='lunar')

#============================================================================#
# Generate total data
# Background = 0
# Constant amplitude and phase
#============================================================================#
data1T = generate_tides('2016-01-01', '2016-03-31', amps=[0,10,10], phase='C',
                        longIncr=15, nRange=[2], sRange=[2], 
                        filename='TT_0bg_CACP.txt', component='s+l')



# Bin by solar local time
means1T = bin_by_solar(data1T, '0bg_CACP')
print('Done with binning by solar')

# Subtract the averages according to solar local time
nosol1T = remove_solar(data1T, '0bg_CACP_slt_bin.txt')
print('Done with removing solar')


compare_with_plot(nosol1T[:,6], data1L[:,6])

# Bin the results by lunar local time
#result = bin_by_lunar(nosol1T, '1T')
#print('Done with binning by lunar')
# The following plot is done WITHOUT binning by LLT just to see how it looks.
#compare_with_plot(result[:,6], data1L[:,6])
#print('Done')