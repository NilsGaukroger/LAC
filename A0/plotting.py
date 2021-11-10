# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 09:05:29 2021

@author: nilsg
"""

# LAC course
# Assignment #0

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

save = True

#%% Part 1: Spanwise values

# read data from .ind file
columns = ['s', 'A', 'AP', 'PHI0', 'ALPHA0', 'U0', 'FX0', 'FY0', 'M0', 'UX0', 'UY0', 'UZ0', 'Twist', 'X_AC0', 'Y_AC0', 'Z_AC0', 'CL0', 'CD0', 'CM0', 'CLp0', 'CDp0', 'CMp0', 'F0', "F'", 'CL_FS0', "CLFS'", 'V_a', 'V_t', 'Tors.', 'vx', 'vy', 'chord', 'CT', 'CP', 'angle', 'v_1', 'v_2', 'v_3']
df = pd.read_csv('46320_LAC\dtu_10mw\DTU_10MW_rigid_hawc2s_u8000.ind',delim_whitespace=True,skiprows=1,names=columns)

# plots
plots = ['A', 'AP', 'CL0', 'CD0', 'CP', 'CT']
fig, axs = plt.subplots(3,2,sharex=True)
fig.tight_layout(w_pad=2.5)
x = [0,0,1,1,2,2]; y = [0,1,0,1,0,1]
ylabels = ['Ax ind factor [-]', 'Tan ind factor [-]', 'Cl [-]', 'Cd [-]', 'CP [-]', 'CT [-]']
for i in range(len(plots)):
    axs[x[i],y[i]].plot(df['s'],df[plots[i]])
    axs[x[i],y[i]].set(ylabel=ylabels[i])
    axs[x[i],y[i]].grid(True)
axs[2,0].set(xlabel='Spanwise position [m]')
axs[2,1].set(xlabel='Spanwise position [m]')
if save == True:
    plt.savefig('figs/Part1_subplots.pdf',bbox_inches='tight')
plt.show()

#%% Part 2: Multiple TSR values

# Determine require rot. speed [RPM]
TSR = np.arange(6.0, 9.5, 0.5)
R = 88.93420444
Uinf = np.arange(7.997, 8.003, 0.001)
omega = TSR * Uinf / R

N = omega * 60 / 2 / np.pi
N[Uinf==8] = 6.442484

# Read data from .ind files
columns2 = ['V', 'P', 'T', 'Cp', 'Ct', 'Pitch Q', 'Flap M', 'Edge M', 'Pitch', 'Speed', 'Tip x', 'Tip y', 'Tip z', 'J_rot', 'J_DT']
power = pd.read_csv('46320_LAC\dtu_10mw\DTU_10MW_rigid_hawc2s.pwr', delim_whitespace=True, skiprows=1, names=columns2)
data = []
for i in range(len(TSR)):
    data.append(pd.read_csv('46320_LAC\dtu_10mw\DTU_10MW_rigid_hawc2s_u' + str(int(Uinf[i]*1000)) + '.ind', delim_whitespace=True, skiprows=1, names=columns))

# plots
colours = plt.cm.viridis(np.linspace(0,0.9,len(TSR)))

plots2 = ['P','Cp','T','Ct']
ylabels2 = ['P [kW]', 'Cp [-]', 'T [kN]', 'Ct [-]']
leg = ['Cp']

fig,axs = plt.subplots(2,1,sharex=True)
ax2 = []
x2 = [0,1]
for i in range(int(len(plots2)/2)):
    axs[x2[i]].plot(TSR,power[plots2[2*i]])
    ax2.append(axs[x2[i]].twinx())
    ax2[i].plot(TSR,power[plots2[2*i+1]],color='red')
    axs[x2[i]].set(ylabel=ylabels2[2*i])
    ax2[i].set(ylabel=ylabels2[2*i+1])
    ax2[i].spines['right'].set_color('red')
    ax2[i].yaxis.label.set_color('red')
    ax2[i].tick_params(axis='y', colors='red')
axs[1].set(xlabel = 'TSR [-]')
if save == True:
    plt.savefig('fig/Part2a_subplots.pdf')
plt.show()

fig, axs = plt.subplots(3,2,sharex=True)
fig.tight_layout(w_pad=2.5)
for i in range(len(plots)):
    for j in range(len(TSR)):
        axs[x[i],y[i]].plot(data[j]['s'],data[j][plots[i]],color=colours[j])
    axs[x[i],y[i]].set(ylabel=ylabels[i])
    axs[x[i],y[i]].grid(True)
axs[2,0].set(xlabel='Spanwise position [m]')
axs[2,1].set(xlabel='Spanwise position [m]')
if save == True:
    plt.savefig('fig/Part2b_subplots.pdf',bbox_inches='tight')
plt.show()

#%% Part 3: Chord, twist, thickness

columns3 = ['r', 'c', 't', 'set', 'opt']
blade = pd.read_csv('46320_LAC\dtu_10mw\data\DTU_10MW_RWT_ae.dat',delim_whitespace=True,skiprows=2,names=columns3)

fig,axs = plt.subplots(3,1,sharex=True)
plots3 = ['c','t']
ylabels3 = ['Chord [m]', 'Thickness [%]', 'Twist [deg]']
for i in range(2):
    axs[i].plot(blade['r'],blade[plots3[i]])
    axs[i].set(ylabel=ylabels3[i])
    axs[i].grid(True)
axs[2].plot(data[4]['s'],np.rad2deg(data[4]['Twist']))
axs[2].set(ylabel=ylabels3[2])
axs[2].set(xlabel='Spanwise position [m]')
axs[2].grid(True)
if save == True:
    plt.savefig('fig/Part3_subplots.pdf',bbox_inches='tight')
plt.show()

#%% Bonus level: Investigate flexibility

# Read data from .ind files
power2 = pd.read_csv('46320_LAC\dtu_10mw\DTU_10MW_flexible_hawc2s.pwr', delim_whitespace=True, skiprows=1, names=columns2)
data2 = []
for i in range(len(TSR)):
    data2.append(pd.read_csv('46320_LAC\dtu_10mw\DTU_10MW_flexible_hawc2s_u' + str(int(Uinf[i]*1000)) + '.ind', delim_whitespace=True, skiprows=1, names=columns))
    
# plots
    
fig, axs = plt.subplots(3,2,sharex=True)
fig.tight_layout(w_pad=2.5)
for i in range(len(plots)):
    axs[x[i],y[i]].plot(data[4]['s'],data[4][plots[i]])
    axs[x[i],y[i]].plot(data2[4]['s'],data2[4][plots[i]])
    axs[x[i],y[i]].set(ylabel=ylabels[i])
    axs[x[i],y[i]].grid(True)
axs[1,1].legend(['Rigid','Flexible'])
axs[2,0].set(xlabel='Spanwise position [m]')
axs[2,1].set(xlabel='Spanwise position [m]')
if save == True:
    plt.savefig('fig/Bonus_subplots.pdf',bbox_inches='tight')
plt.show()

plots3 = ['Cp','Ct']
fig,axs = plt.subplots(2,1,sharex=True)
for i in range(2):
    axs[i].plot(TSR,power[plots3[i]])
    axs[i].plot(TSR,power2[plots3[i]])
    axs[i].grid(True)
    axs[i].set(ylabel=plots3[i]+' [-]')
axs[1].set(xlabel = 'TSR [-]')
axs[0].legend(['Rigid','Flexible'])
if save == True:
    plt.savefig('fig/Bonusa_subplots.pdf')
plt.show()