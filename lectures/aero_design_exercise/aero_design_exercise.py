# -*- coding: utf-8 -*-
"""Example solving the three equations for pointwise blade design
"""
import numpy as np
from scipy.optimize import least_squares

#%% PART 1: Residual

R = 35  # length of blade [m]
tsr = 9.0  # TSR [-]
B = 3  # number of blades
a = 1/3  # axial induction [-]
r = 8  # blade span [m]
t = 1.5  # absolute thickness at blade span [m]
cl = 0.6  # design lift coefficient [-]
clcd = 120  # design cl/cd [-]
alpha = 4  # design AoA [deg]


def part1_residuals(x,res=True):
    """Calculate the residuals on the equations
        Set res=False to return the intermediate variables instead of residuals
    """
    # unpack x
    c, ap = x
    
    # calculate intermediate variables
    phi   = np.arctan(((1-a)/(1+ap)) * (R/(r*tsr)))
    cd    = cl/clcd
    cy    = cl*np.cos(phi) + cd*np.sin(phi)
    cx    = cl*np.sin(phi) - cd*np.cos(phi)
    f     = (B/2) * ((R-r)/(r*np.sin(phi)))
    F     = (2/np.pi) * np.arccos(np.exp(-f))
    sigma = (c*B)/(2*np.pi*r)
    # a     = 1 / (((4*F*np.sin(phi)**2) / (sigma * cy)) + 1)
    CPloc = ((1-a)**2 + (tsr*(r/R)) * (1+ap)**2) * tsr * (r/R) * sigma * cx
    CTloc = ((1-a)**2 + (tsr*(r/R)) * (1+ap)**2) * sigma * cy
    
    # calculate residuals
    res_c = 4*np.pi*r*np.sin(phi)**2*F*((2*a)/(cy*B*(1-a))) - c # residual on chord [m]
    res_ap = 1 / ((4*F*np.sin(phi)*np.cos(phi)/(sigma*cx)) - 1) - ap # residual on tangential induction [-]
    out = [res_c, res_ap]  # return a single output vector
    if res:
        return out
    else:
        return that, cl, cd, phi, sigma, CPloc, CTloc

# test the values
x = [3, 0.001]  # [c (m), a']
out = part1_residuals(x)

print(out)

#%% PART 2: One point with fixed values

x0 = [3, 0.001]  # initial guess [c (m), a']
bounds = ((0, 0), (np.inf, 1))  # lower, upper bounds on solver
c, ap = least_squares(part1_residuals, x0, bounds=bounds).x  # solve the eqns

print(f'c = {c:.4f} m')
print(f'ap = {ap:.4f}')

#%% PART 3: Thickness and design polynomials

def thickness(r):
    """Absolute thickness [m] as a function of radius [m] for 35-m blade"""
    if r < 5:
        r = 5
    t = 9.35996E-08*r**6 - 1.2911E-05*r**5 + 7.15038E-04*r**4 - 2.03735E-02*r**3 + 3.17726E-01*r**2 - 2.65357E+00*r + 1.02616E+01
    return t


def alpha_des(that):
    """Design AOA [deg] as a function of t/c [%]"""
    if that < 36:
    # if 15 <= that < 36:
        alpha_des = 2.32706E-05*that**5 - 2.87331E-03*that**4 + 1.36343E-01*that**3 - 3.10470E+00*that**2 + 3.38460E+01*that - 1.36500E+02
    elif that >= 36:
    # elif 36 <= that <= 100:
        alpha_des = -0.0078*that + 0.7813
    # else:
    #     raise Exception('Could not compute AOA, t/c outside of range 15 to 100%')
    return alpha_des


def cl_des(that):
    """Design cl [-] as a function of t/c [%]"""
    if that < 36:
    # if 15 <= that < 36:
        cl_des = -7.34862E-07*that**5 + 1.10229E-04*that**4 - 6.40432E-03*that**3 + 1.79563E-01*that**2 - 2.43397E+00*that + 1.36000E+01
    elif that >= 36:
    # elif 36 <= that <= 100:
        cl_des = -0.0094*that + 0.9375
    # else:
    #     raise Exception('Could not compute cl, t/c outside of range 15 to 100%')
    return cl_des


def clcd_des(that):
    """Design cl/cd [-] as a function of t/c [%]"""
    if that < 36:
    # if 15 <= that < 36:
        clcd_des = -8.10212E-03*that**4 + 8.73313E-01*that**3 - 3.41180E+01*that**2 + 5.66297E+02*that - 3.24932E+03
    elif that >= 36:
    # elif 36 <= that <= 100:
        clcd_des = -0.8906*that + 89.063
    # else:
    #     raise Exception('Could not compute clcd, t/c outside of range 15 to 100%')
    return clcd_des


import matplotlib.pyplot as plt


fig, axs = plt.subplots(2, 2, num=1, clear=True, figsize=(8, 3.5))

r    = np.linspace(0, 35, 201)
that = np.linspace(15, 100, 201)

# Loop for vectors
thickness_vec = []
cl_des_vec = []
clcd_des_vec = []
alpha_des_vec = []
for i in range(201):
    thickness_vec.append(thickness(r[i]))
    cl_des_vec.append(cl_des(that[i]))
    clcd_des_vec.append(clcd_des(that[i]))
    alpha_des_vec.append(alpha_des(that[i]))

# absolute thickness vs. blade length
ax = axs[0, 0]
ax.plot(r, thickness_vec)
ax.set_xlabel('Blade length [m]'); ax.set_ylabel('Thickness [m]')
ax.grid('on')

# cl,des versus relative thickness
ax = axs[0, 1]
ax.plot(that, cl_des_vec)
ax.set_xlabel('Relative thickness [%]'); ax.set_ylabel('$c_{l,des}$ [-]')
ax.grid('on')

# clcd,des versus relative thickness
ax = axs[1, 0]
ax.plot(that, clcd_des_vec)
ax.set_xlabel('Relative thickness [%]'); ax.set_ylabel('$c_l/c_{d,des}$ [-]')
ax.grid('on')

# alpha,des versus relative thickness
ax = axs[1, 1]
ax.plot(that, alpha_des_vec)
ax.set_xlabel('Relative thickness [%]'); ax.set_ylabel('$\\alpha_{des}$ [deg]')
ax.grid('on')

plt.tight_layout()

#%% PART 4: One span, design polynomials

r = 8

def residuals(x, res=True, tsr=8):
    """Calculate the residuals on the equations. Optionally return BEM values for debugging.
    """
    c, ap = x
    # calculate that, cl and cl/cd from design polynomials
    that = (thickness(r)/c)*100
    cl = cl_des(that)
    clcd = clcd_des(that)
    # intermediate variables
    phi   = np.arctan(((1-a)/(1+ap)) * (R/(r*tsr)))
    cd    = cl/clcd
    cy    = cl*np.cos(phi) + cd*np.sin(phi)
    cx    = cl*np.sin(phi) - cd*np.cos(phi)
    f     = (B/2) * ((R-r)/(r*np.sin(phi)))
    F     = (2/np.pi) * np.arccos(np.exp(-f))
    sigma = (c*B)/(2*np.pi*r)
    # residuals
    res_c = 4*np.pi*r*np.sin(phi)**2*F*((2*a)/(cy*B*(1-a))) - c # residual on chord [m]
    res_ap = 1 / ((4*F*np.sin(phi)*np.cos(phi)/(sigma*cx)) - 1) - ap # residual on tangential induction
    out = [res_c, res_ap]  # return a single output vector
    if res:
        return out
    else:
        return that, cl, cd, phi, sigma

x0 = [7, 0.001]  # initial guess [c (m), a']
out = residuals(x0)
print(out)

bounds = ((0, 0), (np.inf, 1))  # lower, upper bounds on solver
c, ap = least_squares(residuals, x0, bounds=bounds).x  # solve the eqns

print(f'c = {c:.4f} m')
print(f'ap = {ap:.4f}')

#%% PART 5: Multiple spans, design polynomials

radii = np.arange(5, 35, 0.5)  # spanwise points to analyse
out = np.empty((radii.size, 2))  # initialize output array
t = np.empty((radii.size))
for i, r in enumerate(radii):
    t[i] = thickness(r)
    # bounds = ((thickness(r), 0), (thickness(r)/0.15, 1))  # lower, upper bounds on solver
    res = least_squares(residuals, x0, bounds=bounds)  # solve
    out[i] = res.x  # save the solutions
    x0 = res.x  # update initial guess for next radius

print(out[:5])

#%% PART 6: Calculate chord, twist, and relative thickness

# calculate twist and relative thickness
c, ap = out.T
that  = t / c
alpha = np.empty((radii.size))
for i in range(len(that)):
    alpha[i] = alpha_des(that[i]*100)

phi   = np.arctan(((1-a)/(1+ap)) * (R/(radii*tsr)))
theta = np.rad2deg(phi) - alpha

# plot the chord, twist and thickness
fig, axs = plt.subplots(3, 1, num=2, figsize=(7, 6), clear=True)
for i, (val, label) in enumerate([(c, 'Chord [m]'),
                                  (theta, 'Twist [deg]'),
                                  (that, 'Rel. thickness [%]')]):
    axs[i].plot(radii, val)
    axs[i].grid('on'); axs[i].set_ylabel(label)
    if i == 2:
        axs[i].set_xlabel('Radius [m]')
plt.tight_layout()

#%% Question 8: Compare blade designs for different TSRs and induction values

# Varying TSR

tsr_lst = np.arange(6.5, 10, 0.5)
c = np.empty((radii.size,tsr_lst.size))

for j, tsr in enumerate(tsr_lst):
    for i, r in enumerate(radii):
        t[i] = thickness(r) # calculate absolute thickness at each radial point using design polynomial
        res = least_squares(residuals, x0, bounds=bounds, args=[tsr]) # minimise residuals to find c, ap for each radial point
        out[i] = res.x # assign output as result of least squares
        x0 = res.x # initial guess for next spanwise position = result of previous
        
    c, ap = out.T # unpack results
    that = t / c # calculate relative thickness from thickness and chord distributions
    
    # determine angle of attack from relative thickness and design polynomial
    alpha = np.empty((radii.size))
    for i in range(len(that)):
        alpha[i] = alpha_des(that[i]*100)
        
    # determine flow angle and twist
    phi   = np.arctan(((1-a)/(1+ap)) * (R/(radii*tsr)))
    theta = np.rad2deg(phi) - alpha