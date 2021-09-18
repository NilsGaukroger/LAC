# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 14:24:14 2021

@author: nilsg
"""

V1 = 11.4
I1 = 0.16
R1 = 178.3/2
V1_max = 11.4*(1+2*I1)

R2_guess = 100
I2 = 0.14

eps = 0.1
n = 0

while eps > 1e-5:
    V2 = ((V1**3 * R1**2) / (R2_guess**2))**(1/3)
    R2 = (((V1_max)**2 / (V2*(1+2*I2))**2) * R1**2)**(1/2)
    eps = abs(R2 - R2_guess)
    
    R2_guess = R2
    n += 1

print('New diameter converged to {:0.1f}m in {:d} iterations'.format(2*R2,n))