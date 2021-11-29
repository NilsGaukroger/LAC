# -*- coding: utf-8 -*-
"""Plot a statistic versus wind speed.
"""
import matplotlib.pyplot as plt
import numpy as np
from _loads_utils import load_stats, load_hawc2s

hawc2s_path = '../res/redesign/steady/redesign_cont_HS2.pwr'  # path to .pwr or .opt file
stat_dir1 = '../res/dtu10mw/turb/'  # stats file, turbine 1  !!! END WITH SLASH !!!
stat_dir2 = '../res/redesign/turb/'  # stats file, turbine 2  !!! END WITH SLASH !!!
labels = ['DTU 10 MW, TC A', 'Redesign, TC B']  # legend labels
labels2 = ['Redesign, TC B']

# values to plot and indices of that value for the two turbines
wind_idxs = [15, 15]  # indices of vy wind for turbine 1 and turbine 2
#               Label                  Idx_t1  Idx_t2
plot_vals = [['Pitch angle [deg]',       4,    4],
             ['Rotor speed [rad/s]',     10,   10],
             ['Thrust [kN]',             13,   13],
             ['AoA @ 0.84R m',           None,   63],
             ['Cl @ 0.84R m',            None,   66],
             ['Tower-base FA [kNm]',     19,   19],
             ['Tower-base SS [kNm]',     20,   20],
             ['Yaw-bearing pitch [kNm]', 22,   22],        
             ['Yaw-bearing roll [kNm]',  23,   23],
             ['Shaft torsion [kNm]',     27,   27],
             ['OoP BRM [kNm]',           28,   28],
             ['IP BRM [kNm]',            29,   29],
             ['Generator torque [Nm]',   72,   72],
             ['Electrical power [W]',    102,  102], 
             ['Tower clearance [m]',     110,  110], 
             ]

fig_names = ['pitchAngle','rotorSpeed','thrust','AoA','Cl','TBFA','TBSS','YBP','YBR','ST','OoP_BRM','IP_BRM','genTrq','elec_power','tower_clear']

# load the stuff we need from the HAWC2S .pwr file for the operational data comparisons
h2s_u, h2s_pitch, h2s_rotspd, h2s_paero, h2s_thrust, h2s_aerotrq = load_hawc2s(hawc2s_path)

geneff = 0.895
h2s_paero = h2s_paero * geneff
h2s_aerotrq = h2s_aerotrq * geneff
h2s_thrust = h2s_thrust * 1e-3
hawc2s_vars = [h2s_pitch,h2s_rotspd,h2s_thrust,h2s_aerotrq,h2s_paero]

# hawc2s_idxs = [0,1,2,12,13]
hawc2s_idxs = np.array([4,10,13,72,102])

# =======================================================================================
# you shouldn't need to change below this line :)

axs = []  # initialize list of axes
# colors = ['#1f77b4', '#ff7f0e']  # colors for turbines 1 and 2
colors = ['0.3', 'r']  # colors for turbines 1 and 2
markers = ['x', '.']  # markers for turbines 1 and 2

# load the stats data for turbines 1 and 2 (it's ugly i'm sorry)
_, idxs_1, means_1 = load_stats(stat_dir1 + 'stats_mean.txt')  # means
_, _, mins_1 = load_stats(stat_dir1 + 'stats_min.txt')  # mins
_, _, maxs_1 = load_stats(stat_dir1 + 'stats_max.txt')  # maxs
wind_1 = means_1[:, idxs_1 == wind_idxs[0]]
_, idxs_2, means_2 = load_stats(stat_dir2 + 'stats_mean.txt')  # means
_, _, mins_2 = load_stats(stat_dir2 + 'stats_min.txt')  # mins
_, _, maxs_2 = load_stats(stat_dir2 + 'stats_max.txt')  # maxs
wind_2 = means_2[:, idxs_2 == wind_idxs[1]]

# plotting settings for that turbine
# m = markers[iturb]; mec = colors[iturb]

# loop over the channels to plot
for i, (ylabel, idx_t1, idx_t2) in enumerate(plot_vals):
    
    # make the plot and initialize some plot settings
    fig, ax = plt.subplots(num=1 + i, figsize=(7, 3), clear=True)
    ax.grid('on')
    ax.set(xlabel='Wind speed [m/s]', ylabel=ylabel)

    # loop through the turbines
    handles = []
    for iturb, (idxs, means, mins, maxs, wind) in enumerate(zip(
            [idxs_1, idxs_2], [means_1, means_2], [mins_1, mins_2],
            [maxs_1, maxs_2], [wind_1, wind_2])):
        
        # get the channel idx, skip if None
        chan_idx = [idx_t1, idx_t2][iturb]
        if chan_idx is None:
            continue

        # isolate the channel to plot
        meanval = means[:, idxs == chan_idx]
        minval = mins[:, idxs == chan_idx]
        maxval = maxs[:, idxs == chan_idx]
    
        # get plotting marker and color
        m = markers[iturb]; mec = colors[iturb]
    
        # only plot minimum value if tower clearance
        if 'clearance' in ylabel:
            l, = ax.plot(wind, minval, marker='.', mec=mec, mfc='none', alpha=0.8, linestyle='none')
        # otherwise plot mean, max and min
        else:
            l, = ax.plot(wind, meanval, marker='.', mec=mec, mfc='none', alpha=0.8, linestyle='none')
            ax.plot(wind, minval, marker='.', mec=mec, mfc='none', alpha=0.8, linestyle='none')
            ax.plot(wind, maxval, marker='.', mec=mec, mfc='none', alpha=0.5, linestyle='none')
            # ax.plot(wind, minval, marker='x', mec=mec, ms=4, mfc='none', alpha=0.8, linestyle='none')
            # ax.plot(wind, maxval, marker='x', mec=mec, ms=4, mfc='none', alpha=0.5, linestyle='none')
        handles.append(l)
        if iturb == 1 and (idx_t2 in hawc2s_idxs):
            l2, = ax.plot(h2s_u,hawc2s_vars[np.where(hawc2s_idxs==idx_t2)[0][0]],'--',color='black')
            handles.append(l2)
            labels.append('HAWC2S')
        if iturb == 1 and ('Cl' in ylabel):
            ax.hlines(1.8140, 4, 25, color="black", linestyle='dashed')
        if iturb == 1 and ('AoA' in ylabel):
            ax.hlines(16, 4, 25, color="black", linestyle='dashed',label=r'$\alpha_{max}$')
        if 'Tower-base' in ylabel:
            ax.axvspan(6.5, 9.5, color='orange', alpha=0.3)
    
    if 'AoA' in ylabel or 'Cl' in ylabel:
        ax.legend(labels=labels2)
    else:
        ax.legend(handles=handles, labels=labels)
    fig.tight_layout()
    
    # save figures to /figs/ folder
    plt.savefig('../figs/Part1and2/' + fig_names[i] + '.png')
    plt.show()