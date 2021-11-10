# -*- coding: utf-8 -*-
"""Plot the Campbell diagram for Ronnie
"""
import matplotlib.pyplot as plt
import numpy as np
from _ronnie_utils import get_nr_matrices, get_eig


# define our input variables
Omegas = np.linspace(0, 1, 21)  # rotational frequencies to plot [rad/s]

# =======================================================================================
# you shouldn't need to change anything below this line :)

kwargs = dict(dofs=[1, 1, 1, 1, 0],  # flags to enable degrees of freedom [Bld1, Bld2, Bld3, x, y]
              M=50_000,  # equivalent rotor/nacelle/tower mass [N/m]
              kx=200_000,  # equivalent tower stiffness [N/m]
              m=500,  # equivalent blade mass [kg]
              l=30,  # equivalent blade length [m]
              omega0=2*np.pi*0.8)  # blade natural frequency [rad/s]

# initialize the output
f_nr = np.empty((Omegas.size, sum(kwargs['dofs'])))
for iom, Omega in enumerate(Omegas):
    
    # get the system matrices in the non-rotating frame
    M_nr, C_nr, K_nr = get_nr_matrices(Omega, **kwargs)
    n = M_nr.shape[0]  # number DOFs
    
    # get the eigenvalues and eigenvectors
    [L, P] = get_eig(M_nr, C_nr, K_nr)
    
    # isolate the states that correspond to displacements
    L = L[::2]
    P = P[:n, ::2]
    
    # calculate natural frequencies and damping
    omegas = np.abs(L)  # rad/s
    freqs = omegas / 2 / np.pi  # Hz
    zetas = -np.cos(np.angle(L)) * 100  # % critical
    
    # calculate the symmetric, forwards whirling and backwards wirling magnitudes
    a0, a1, b1 = P[:3, :]  # extract the blade components from eigenvectors
    re, im = np.real, np.imag  # assign aliases to functions for convenience
    S_mag = np.abs(a0)  # symmetric is just a0
    BW_mag = 0.5*np.sqrt((im(a1)-re(b1))**2 + (re(a1)+im(b1))**2)  # BW, eqn in notes
    FW_mag = 0.5*np.sqrt((im(a1)+re(b1))**2 + (re(a1)-im(b1))**2)  # FW, eqn in notes
    x_mag = np.abs(P[3])  # center mass displacement
    
    # sort mode shapes into fore-aft, BW, S, FW
    i_x = x_mag.argmax()
    i_BW = BW_mag.argmax()
    i_S = S_mag.argmax()
    i_FW = FW_mag.argmax()
    if len(set((i_x, i_BW, i_S, i_FW))) < n:
        if Omega > 0:
            print(f'WARNING! Cannot sort modes for Omega = {Omega:.2f}.')
        freqs = np.sort(freqs)
        i_x, i_BW, i_S, i_FW = range(4)
    f_nr[iom] = freqs[[i_x, i_BW, i_S, i_FW]]


# make the campbell diagram
fig = plt.figure(1, figsize=(7, 3.8), clear=True)
plt.plot(Omegas, f_nr, 'o')
plt.plot(Omegas, Omegas/2/np.pi, ':', c='0.2')
plt.plot(Omegas, 3*Omegas/2/np.pi, '--', c='0.2')
plt.plot(Omegas, kwargs['omega0']/2/np.pi*np.ones_like(Omegas), ':', c='0.5')
plt.legend(labels=['Tower side-side', 'Backwards whirling', 'Symmetric', 'Forwards whirling',
                   '1P', '3P'],
           bbox_to_anchor=(0., 1.04, 1., .102), loc='lower left',
           ncol=3, mode='expand', borderaxespad=0., fontsize=10)
plt.grid('on')
plt.xlabel('Rotational frequency [rad/s]')
plt.ylabel('Frequencies [Hz]')
plt.tight_layout()

