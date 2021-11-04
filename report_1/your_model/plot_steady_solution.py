# -*- coding: utf-8 -*-
"""An exercise to plot the mean loads for a series of channels and comapre to the theory.
"""
import matplotlib.pyplot as plt
import numpy as np
from _loads_utils import load_stats, load_hawc2s


hawc2s_path = './dtu10mw_res/DTU_10MW_hawc2s.pwr'  # path to .pwr or .opt file
stats_path = './dtu10mw_res/dtu10mw_steady/stats_mean.txt'  # path to mean steady stats

dz_tb = 119  # distance from hub center to tower base [m]
dz_yb = 3.37  # distance from hub center to yaw bearing [m]
geneff = 0.94  # generator/gearboox efficienty [%]
Mgrav = 5750  # yaw-bearing pitch moment due to gravity [kNm]

# define the channels and names to plot
channels = {4: 'Pitch angle [deg]',
            10: 'Rotor speed [rad/s]',
            13: 'Thrust [kN]',
            72: 'Generator torque [Nm]',
            102: 'Electrical power [W]',
            63: 'Angle of attack @ 2/3 R [deg]',
            66: 'Cl @ 2/3 R [-]',
            19: 'Tower-base FA [kNm]',
            20: 'Tower-base SS [kNm]',
            22: 'Yaw-bearing pitch [kNm]',
            23: 'Yaw-bearing roll [kNm]',
            27: 'Shaft torsion [kNm]',
            28: 'Out-of-plane BRM [kNm]',
            29: 'In-plane BRM [kNm]'}
i_wind = 15  # wind channel, needed for plotting versus wind speed

# load the HAWC2 data from the stats file
files, idxs, data = load_stats(stats_path)
wind = data[:, idxs == i_wind]

# load the stuff we need from the HAWC2S .pwr file for the operational data comparisons
h2s_u, h2s_pitch, h2s_rotspd, h2s_paero, h2s_thrust, h2s_aerotrq = load_hawc2s(hawc2s_path)

# loop over each channels and plot the steady state with the theory line
for iplot, (ichan, name) in enumerate(channels.items()):
    
    # extract hawc2 wind and channel to plot from the HAWC2 stats
    h2_wind = data[:, idxs == i_wind]  # wind [m/s]
    HAWC2val = data[:, idxs == ichan]

    # hawc2 channels we need for the theoretical calculations
    h2_thrust = data[:, idxs == 13]  # thrust [kN]
    h2_aero_trq = data[:, idxs == 72] / geneff * 1e-3  # aerodynamic torque [kNm]
    
    # PART 1. Theoretical lines. Taken from HAWC2S pwr file!
    if ichan == 4:  # pitch angle
        u_theory = h2s_u
        theory = h2s_pitch
    elif ichan == 10:  # rotor speed
        u_theory = h2s_u
        theory = h2s_rotspd
    elif ichan == 13:  # thrust
        u_theory = h2s_u
        theory = h2s_thrust * 1e-3
    elif ichan == 72:  # generator torque
        u_theory = h2s_u
        theory = h2s_aerotrq * geneff
    elif ichan == 102:  # electrical power
        u_theory = h2s_u
        theory = h2s_paero * geneff

    # aerodynamic parameters
    elif ichan == 63:  # angle of attack
        u_theory = h2_wind
        theory = np.nan * np.ones_like(u_theory)
    elif ichan == 66:  # lift coefficient
        u_theory = h2_wind
        theory = np.nan * np.ones_like(u_theory)

    # PART 2. Theoretical lines. Equations in lecture, calculated with hawc2 channels.
    elif ichan == 19:  # tower-base fore-aft
        u_theory = h2_wind
        theory = h2_thrust * dz_tb - Mgrav
    elif ichan == 20:  # tower-base side-side
        u_theory = h2_wind
        theory = h2_aero_trq
    elif ichan == 22:  # yaw bearing pitch
        u_theory = h2_wind
        theory = h2_thrust * dz_yb - Mgrav
    elif ichan == 23:  # yaw bearing roll
        u_theory = h2_wind
        theory = h2_aero_trq
    elif ichan == 27:  # shaft torsion
        u_theory = h2_wind
        theory = -h2_aero_trq
    elif ichan == 28:  # blade root out-of-plane moment
        u_theory = h2_wind
        theory = h2_thrust * -60 / 3
    elif ichan == 29:  # blade root in-plane moment
        u_theory = h2_wind
        theory = h2_aero_trq/ 3

    # plot the results
    fig = plt.figure(1 + iplot, figsize=(7, 3), clear=True)
    plt.plot(u_theory, theory, '--', c='0.2')  # theoretical line
    plt.plot(wind, HAWC2val, 'o')  # HAWC2 steady results
    plt.grid('on')
    plt.xlabel('Wind speed [m/s]')
    plt.ylabel(name)
    plt.tight_layout()
