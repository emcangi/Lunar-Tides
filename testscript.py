# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 17:48:40 2016

@author: emc
"""

from solar_extraction import *

minus_sol = []
totals = []
results = []
dts = [0.25, 0.5, 1]
binsize = 1

for dt in dts:
    # ==========================================================================
    # Generate SW2 data
    # Background = 0
    # Constant amplitude and phase
    # ==========================================================================
    data1S, SW2_120 = generate_tides('2016-01-01', '2016-01-15',
                                     amps=[0, 10, 10], phase='C', dt=dt,
                                     nrange=[2], srange=[2],
                                     filename='SW2_revised_dt{}.txt'.format(dt),
                                     component='solar')

    # ==========================================================================
    # Generate M2 data
    # Background = 0
    # Constant amplitude and phase
    # ==========================================================================
    data1M, M2_120 = generate_tides('2016-01-01', '2016-01-15',
                                    amps=[0,10,10], phase='C', dt=dt,
                                    nrange=[2], srange=[2],
                                    filename='M2_revised_dt{}.txt'.format(dt),
                                    component='lunar')

    # ==========================================================================
    # Generate total data
    # Background = 0
    # Constant amplitude and phase
    # ==========================================================================
    data1T, TT_120 = generate_tides('2016-01-01', '2016-01-15',
                                    amps=[0,10,10], phase='C', dt=dt,
                                    nrange=[2], srange=[2],
                                    filename='TT_revised_dt{}.txt'.format(dt),
                                    component='s+l')

    # Save a file for one longitude only
    # np.savetxt('SW2_revised_lon=-120_dt{}.txt'.format(dt), SW2_120,
    #            fmt='%20.4f', delimiter='\t')
    # np.savetxt('M2_revised_lon=-120.txt', M2_120, fmt='%20.4f',
    #            delimiter='\t')
    # np.savetxt('TT_revised_lon=-120.txt', TT_120, fmt='%20.4f',
    #            delimiter='\t')

    # BIN GENERATED DATA BY SOLAR LOCAL TIME ===================================
    means = bin_by_solar(data1T, binsize)
    # np.savetxt('means_dt={}.txt'.format(dt), means, fmt='%-20.4f',
    #             delimiter='\t')


    # SUBTRACT AVERAGES BY SOLAR LOCAL TIME ====================================

    avgs, nosol = remove_solar(data1T, means, binsize)
    # longs = np.around(avgs[:, 2], decimals=4)
    # subs120 = avgs[np.where(longs==-120)]
    # np.savetxt('SW2_lon=-120_dt={}_subtracts_restricted_to_120.txt'.format(dt),
    #            avgs, fmt='%-20.4f', delimiter='\t')
    #
    #subs.append(avgs)


    # Bin the results by lunar local time

    means_llt = bin_by_lunar(nosol, binsize)
    final_result = insert_llt_avgs(nosol, means_llt, binsize)

    # build lists of result arrays to use for plotting
    minus_sol.append(nosol)
    totals.append(data1M)
    results.append(final_result)

# plot_vs_date_multi(minus_sol, -120, title='Total SW2+M2 and total - SLT avg, '
#                                       'half lunar cycle, binsz=30 min',
#                    dts=dts, data2=totals, c=['blue', 'deepskyblue'],
#                    lb=['Total - solar avg', 'Original M2'])



plot_vs_date_multi(results, -120, title='Original and reconstructed M2 after '
                                        'LLT binsz, half lunar cycle, binsz=1 hr',
                                        dts=dts, data2=totals,
                                        c=['blue', 'deepskyblue'],
                                        lb=['Reconstructed M2', 'Original M2'])