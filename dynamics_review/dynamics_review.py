# -*- coding: utf-8 -*-
"""A review of modal analysis for linear, time-invariant systems.
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eig


# system parameters (mimicking DTU 10 MW)
m1 = 551_556  # (nacelle mass + hub mass) [kg]
m2 = 122_442  # (rotor mass - hub mass) [kg]
c1 = 6_866  # equivalent tower damping [kg*s]
c2 = 3_599  # equivalent tower damping [kg*s]
k1 = 1_745_883  # equivalent tower stiffness [N/m]
k2 = 1_503_050  # equivalent rotor stiffness [N/m]

# assemble the mass, stiffness and damping matrices
M = np.array([[m1, 0], [0, m2]])
C = np.array([[c1+c2, -c2], [-c2, c2]])
K = np.array([[k1+k2, -k2], [-k2, k2]])

# assemble the state-space matrices
I, O = np.eye(2), np.zeros((2, 2))
A = np.block([[O, I],
              [-K, -C]])
B = np.block([[I, O],
              [O, M]])


def plot_modes(omegas, V, t=np.linspace(0, 4, num=201), num=1):
    """Plot the oscillations of the two masses for the two mode shapes.
    """
    fig, axs = plt.subplots(1, 2, num=num, figsize=(7, 2.5), clear=True)
    # plot the oscillating undamped mode
    for i in range(2):  # loop over modes
        ax = axs[i]  # select the axes of interest
        mode = np.real(V[:, i] * np.exp(1j * omegas[i] * t)[:, None])  # oscillation versus time
        ax.plot(t, mode)
        ax.grid('on')
        ax.set_title(f'Mode {i+1}', fontsize=10)
        ax.set_xlim([0, t.max()])
        ax.set_xlabel('Time [s]')
    axs[0].legend(labels=['$x_1$', '$x_2$'])
    fig.tight_layout()
    return fig, axs


#%% PART 1. Solving the undamped system (no state space required).

# solve generalized eigenvalue problem
[D,V] = eig(-K,M)

# sort eigenvalues and vectors in increasing order of frequency
argsort = np.argsort(np.abs(D))
D = D[argsort]
V = V[:, argsort]

# calculate frequencies
omegas = abs(np.sqrt(D))    # Omega is n rad/s
freqs  = omegas / (2*np.pi) # freqs in Hz

print('\nUNDAMPED\n')
print('    Mode 1    Mode 2')
print(f'fk  {freqs[0]:.2f}   {freqs[1]:.2f}')
print(f'x1  {V[0, 0]:.2f}   {V[0, 1]:.2f}')
print(f'x2  {V[1, 0]:.2f}   {V[1, 1]:.2f}')

# plot mode oscillation vs. time
fig, axs = plot_modes(omegas, V, num=1)
fig.suptitle('Undamped')


#%% PART 2. Solving the damped system using state space.

# solve the generalized eigenvalue problem
[L, P] = eig(A, B)

# isolate the states that correspond to displacements
L = L[::2]
P = P[:2, ::2]

# sort in increasing order of frequency
argsort = np.argsort(np.abs(L))
L = L[argsort]
P = P[:, argsort]

# calculate natural frequencies and damping
omegas = abs(L)                   # Omega is n rad/s
freqs  = omegas / (2*np.pi)       # freqs in Hz
zetas  = -np.cos(np.angle(L))*100 # If you want to convert zeta to percent critical, multiply by 100!

# intermediates: magnitude and phase of mode
mags = np.abs(P)
pha = np.angle(P)

print('\nDAMPED\n')
print('    Mode 1    Mode 2')
print(f'fk  {freqs[0]:.2f}   {freqs[1]:.2f}')
print(f'zk  {zetas[0]:.2f}   {zetas[1]:.2f}')
print(f'x1  {mags[0, 0]:.2f}<{pha[0, 0]:.2f}   {mags[0, 1]:.2f}<{pha[0, 1]:.2f}')
print(f'x2  {mags[1, 0]:.2f}<{pha[1, 0]:.2f}   {mags[1, 1]:.2f}<{pha[1, 1]:.2f}')

# plot mode oscillation vs. time
fig, axs = plot_modes(omegas, P, num=2)
fig.suptitle('Damped')

