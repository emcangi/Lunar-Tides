# -*- coding: utf-8 -*-
"""
A simple, small script to generate separate solar and lunar tides and then
show them overlapped on a plot for demonstration of similar periods.

Author: Eryn Cangi
LASP REU, CU Boulder
13 June - 29 July, 2016
"""

from lunar_tide_extraction import *

solar = generate_tides('2016-06-01', '2016-06-30', amps=[0, 75, 75],
                        component='solar')

lunar = generate_tides('2016-06-01', '2016-06-30', amps=[0, 75, 75],
                       component='lunar')

plot_vs_date(solar, -105, 'Solar and lunar tides', data2=lunar,
             c=['orange', 'deepskyblue'], lb=['Solar Tide', 'Lunar tide'],
             mode='both')
