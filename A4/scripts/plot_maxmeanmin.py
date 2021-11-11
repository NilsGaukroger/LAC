# -*- coding: utf-8 -*-
"""Plot a statistic versus wind speed.
"""
import matplotlib.pyplot as plt
from _loads_utils import load_stats


stat_dir = '../res/redesign/turb/'  # results directory with statistics files  !!! END WITH SLASH !!!
i_wind = 15  # channel number with the wind speed
i_plot = [4, 10, 13, 15, 19, 20, 22, 23, 27, 28, 29, 72, 102, 110]  # channel indices in .sel file that you want to plot

# dictionary to map .sel index to ylabel for the plot
ylabels = {4: 'Pitch angle [deg]',
           10: 'Rotor speed [rad/s]',
           13: 'Thrust [kN]',
           15: 'Wind speed [m/s]',
           19: 'Tower-base FA [kNm]',
           20: 'Tower-base SS [kNm]',
           22: 'Yaw-bearing pitch [kNm]',
           23: 'Yaw-bearing roll [kNm]',
           27: 'Shaft torsion [kNm]',
           28: 'OoP BRM [kNm]',
           29: 'IP BRM [kNm]',
           72: 'Generator torque [Nm]',
           102: 'Electrical power [W]',
           110: 'Tower clearance [m]'}

# =======================================================================================
# you shouldn't need to change below this line :)

# load the statistics data
_, idxs_means, data_means = load_stats(stat_dir + 'stats_mean.txt')  # means
_, idxs_mins, data_mins = load_stats(stat_dir + 'stats_min.txt')  # mins
_, idxs_maxs, data_maxs = load_stats(stat_dir + 'stats_max.txt')  # maxs
wind = data_means[:, idxs_means == i_wind]


for i, chan_idx in enumerate(i_plot):

    # define ylabel
    ylabel = ylabels[chan_idx]

    # isolate the channel to plot
    meanval = data_means[:, idxs_means == chan_idx]
    minval = data_mins[:, idxs_mins == chan_idx]
    maxval = data_maxs[:, idxs_maxs == chan_idx]
    
    # make the plot
    fig, ax = plt.subplots(num=1 + i, figsize=(7, 3), clear=True)
    # only plot minimum value if tower clearance
    if 'clearance' in ylabel:
        ax.plot(wind, minval, 'o', mfc='w', alpha=0.8)
    # otherwise plot mean, max and min
    else:
        l, = ax.plot(wind, meanval, 'o', mfc='w', alpha=0.8)
        ax.plot(wind, minval, 'x', c=l.get_color(), alpha=0.8)
        ax.plot(wind, maxval, 'x', c=l.get_color(), alpha=0.8)
    ax.grid('on')
    ax.set(xlabel='Wind speed [m/s]', ylabel=ylabel)
    plt.tight_layout()

plt.show()
