# -*- coding: utf-8 -*-
"""Calculate extreme loads.

Use method (a) in IEC 61400-1 (2019), which is extreme times 1.35. Note that we are
incorrectly applying this to all loads, when it is only technically allowed for blade
root moments and deflections.
"""
import re
import matplotlib.pyplot as plt
import numpy as np
from _loads_utils import load_stats



stat_dir = '../res/redesign/turb/'  # results directory with statistics files  !!! END WITH SLASH !!!
i_plot = [19, 20, 22, 23, 27, 28, 29]  # channel indices in .sel file that you want to process
i_wind = 15  # channel number with the wind speed
# undef_tclear = 18.25  # undeflected-blade tower clearance [m]
v_ref = 50  # reference wind speed based on wind class (I=50, 2=42.5, 3=37.5)


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

# figure file names - moved to plot_dels_2turbines
# fig_names = ['TBFA','TBSS','YBP','YBR','ST','OoP_BRM','IP_BRM']

# load the del 4 statistics
stat_file = stat_dir + 'stats_del4.txt'
files, idxs, data_4 = load_stats(stat_file)

# load the del 10 statistics
stat_file = stat_dir + 'stats_del10.txt'
files, idxs, data_10 = load_stats(stat_file)

# get the mean wind for plotting
stat_file = stat_dir + 'stats_mean.txt'
files, idxs, data_mean = load_stats(stat_file)
wind = data_mean[:, idxs == i_wind].squeeze()

# extract the set wind speed value from the filename using regex tricks
wsps = [float(re.findall('[0-9]{1,2}[.][0-9]', f)[0]) for f in files]

# Preallocate list for saving DEL lifetime values
del_life_lst = np.array([])

# loop through the channels
for i, chan_idx in enumerate(i_plot):

    # define ylabel
    ylabel = ylabels[chan_idx]

    # determine which DEL to select
    if 'BRM' in ylabel:
        data = data_10[:, idxs == chan_idx].squeeze()
        m = 10  # 10 exponent for composites
    else:
        data = data_4[:, idxs == chan_idx].squeeze()
        m = 3  # 4 exponent for metals

    # combine short-term dels in given wind speed bin to single value for that bin
    wsp_uniqe = np.unique(wsps)
    st_dels = np.empty(wsp_uniqe.size)
    for j, vj in enumerate(wsp_uniqe):
        wsp_dels = data[np.isclose(wsps, vj)]  # short-term dels
        p = 1 / wsp_dels.size  # probability of each wsp in this bin
        st_dels[j] = sum(p * wsp_dels**m)**(1/m)

    # plot short-term dels versus wind speed
    fig = plt.figure(1 + i, figsize=(7, 3), clear=True)
    plt.plot(wind, data, 'o')
    plt.grid('on')
    plt.xlabel('Wind speed [m/s]')
    plt.ylabel(ylabel)
    plt.tight_layout()

    # for fun, plot the wind-speed-averaged DELs on top
    plt.plot(wsp_uniqe, st_dels, 'or', mec='0.2', ms=7, alpha=0.9)
    
    # save figures as pdf to figs/ folder - moved to plot_dels_2turbines
    # plt.savefig('../figs/Part3/' + fig_names[i] + '.pdf')

    # calculate the lifetime damage equivalent load
    v_ave = 0.2 * v_ref  # average wind speed per IEC 61400-1
    dvj = wsp_uniqe[1] - wsp_uniqe[0]  # bin width. assumes even bins!
    probs = (np.exp(-np.pi*((wsp_uniqe - dvj/2) / (2*v_ave))**2)
             - np.exp(-np.pi*((wsp_uniqe + dvj/2) / (2*v_ave))**2))  # prob of wind in each bin
    del_life = sum(probs * st_dels**m)**(1/m)  # sum DELs ignoring reference number of cycles
    del_life_lst = np.append(del_life_lst,del_life)

    print(ylabel, f'{del_life:.6e}')

#%% Polar plot for comparison with DTU 10 MW
dtu_del_life_lst = np.array([47792.2, 22188.9, 11580.7, 1491.17, 1923.16, 18714.4, 17015.6])
redesign_del_life_lst_normalised = del_life_lst / dtu_del_life_lst
redesign_del_life_lst_normalised = np.append(redesign_del_life_lst_normalised,None)

theta = np.arange(0,2*np.pi,(2*np.pi)/len(redesign_del_life_lst_normalised))
theta = np.append(theta,2*np.pi)
redesign_del_life_lst_normalised = np.append(redesign_del_life_lst_normalised,redesign_del_life_lst_normalised[0])
circle_theta = np.linspace(0,2*np.pi,100)
circle_r     = np.ones((100))

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(circle_theta,circle_r,'--',color='gray',linewidth=2)
ax.plot(theta,redesign_del_life_lst_normalised,'o-',linewidth=2)
ax.set_xticklabels(["TbFA","TbSS","YbTI","YbRI","Shft","FlBR","EdBR","TClear"])
# ax.set_rmin(0.2)
ax.set_rmax(2.5)
ax.set_rticks(np.arange(0.5,3.0,0.5))
ax.set_title('Fatigue loads')
plt.savefig('../figs/Part4/polar_fat.png')
plt.show()