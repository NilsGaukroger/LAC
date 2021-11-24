# -*- coding: utf-8 -*-
"""Calculate extreme loads.

Use method (a) in IEC 61400-1 (2019), which is mean extreme times 1.25*1.35. Note that
we are incorrectly applying this to all loads, when it is only technically allowed for
blade root moments and deflections.
"""
import re
import matplotlib.pyplot as plt
import numpy as np
from _loads_utils import load_stats


stat_dir = '../res/redesign/turb/'  # results directory with statistics files  !!! END WITH SLASH !!!
i_plot = [19, 20, 22, 23, 27, 28, 29, 110]  # channel indices in .sel file that you want to process
i_wind = 15  # channel number with the wind speed

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

# load the min statistics
stat_file = stat_dir + 'stats_min.txt'
files, idxs, data_min = load_stats(stat_file)

# load the mean statistics
stat_file = stat_dir + 'stats_mean.txt'
files, idxs, data_mean = load_stats(stat_file)
wind = data_mean[:, idxs == i_wind].squeeze()

# load the max statistics
stat_file = stat_dir + 'stats_max.txt'
files, idxs, data_max = load_stats(stat_file)

# extract the set wind speed value from the filename using regex tricks
wsps = [float(re.findall('[0-9]{1,2}[.][0-9]', f)[0]) for f in files]

extr_design_load_lst = []

# loop through the channels
for i, chan_idx in enumerate(i_plot):

    # define ylabel
    ylabel = ylabels[chan_idx]

    # isolate the channels to plot
    minval = data_min[:, idxs == chan_idx]
    maxval = data_max[:, idxs == chan_idx]

    # get mean of the extremes for each wind speed bin
    wsp_unique = np.unique(wsps)
    mean_min = np.empty(wsp_unique.size)
    mean_max = np.empty(wsp_unique.size)
    for j, vj in enumerate(wsp_unique):
        mean_min[j] = minval[np.isclose(wsps, vj)].mean()
        mean_max[j] = maxval[np.isclose(wsps, vj)].mean()
    
    # calculate design extreme loads
    if 'clearance' in ylabel.lower():  # tower clearance is weird: no safety factors and min of min
        maxval = np.full_like(maxval, np.nan)  # set maxes to nans, because those are irrelevant
        meanval = np.full_like(maxval, np.nan)  # set maxes to nans, because those are irrelevant
        mean_max[:] = np.nan
        extr_design_load = np.min(minval)  # transform back to tower clearance
        extr_design_load_lst.append(extr_design_load)
    else:
        extremes = np.hstack((mean_min, mean_max)).squeeze()  # combine min and max extremes
        i_ext = np.nanargmax(np.abs(extremes))  # find index of value with max abs value
        fc = 1.35 * extremes[i_ext]  # characteristic load
        extr_design_load = 1.25 * fc  # extreme design load
        extr_design_load_lst.append(extr_design_load)
    print(ylabel, f'{extr_design_load:.6e}')

    # make the plot
    fig = plt.figure(1 + i, figsize=(7, 3), clear=True)
    plt.plot(wind, minval, 'o')  # minimum values
    plt.plot(wsp_unique, mean_min, 'or', mec='0.2', ms=7, alpha=0.8)  # plot mean extremes for fun
    plt.plot(wind, maxval, 'o')  # maximum values
    plt.plot(wsp_unique, mean_max, 'or', mec='0.2', ms=7, alpha=0.8)  # plot mean extremes for fun
    plt.plot(wind, extr_design_load * np.ones_like(wind), lw=2, c='0.2')  # extreme design load
    plt.grid('on')
    plt.xlabel('Wind speed [m/s]')
    plt.ylabel(ylabel)
    plt.tight_layout()

plt.show()

#%% Polar plot for comparison with DTU 10 MW
dtu_extr_design_load_lst = np.array([370603.40625, 110181.290625, 52446.0375, 24233.090625, -23052.09375, -62707.921875, 33696.534375, 6.59303])
redesign_extr_design_load_normalised = extr_design_load_lst / dtu_extr_design_load_lst

theta = np.arange(0,2*np.pi,(2*np.pi)/len(redesign_extr_design_load_normalised))
theta = np.append(theta,2*np.pi)
redesign_extr_design_load_normalised = np.append(redesign_extr_design_load_normalised,redesign_extr_design_load_normalised[0])
circle_theta = np.linspace(0,2*np.pi,100)
circle_r     = np.ones((100))

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(circle_theta,circle_r,'--',color='gray',linewidth=2)
ax.plot(theta,redesign_extr_design_load_normalised,'o-',linewidth=2)
ax.set_xticklabels(["TbFA","TbSS","YbTI","YbRI","Shft","FlBR","EdBR","TClear"])
ax.set_rmin(0.2)
ax.set_rmax(2)
ax.set_rticks(np.arange(0.5,2.5,0.5))
ax.set_title('Extreme loads')
plt.savefig('../figs/Part4/polar_extr.png')
plt.show()