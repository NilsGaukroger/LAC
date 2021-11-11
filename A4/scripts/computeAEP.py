# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 13:48:55 2021

@author: nilsg
"""

import numpy as np

def computeAEP(Umean,Ustd,reliability,pwrs,wsps):
    C = (2/np.sqrt(np.pi))*Umean
    k = (Ustd/Umean)**(-1.086)
    prob = []
    pP   = []
    for i, wsp in enumerate(wsps):
        prob.append(np.exp(-(wsp[0]/C)**k) - np.exp(-(wsp[1]/C)**k))
        pP.append(prob[i]*pwrs[i])
        
    AEP = np.sum(pP) * reliability * 365.25 * 24
    
    return AEP

Umean = 8
Ustd  = 5
reliability = 0.95
pwrs = [0,20,35,40,45,0]
wsps = [[0,5],[5,8],[8,12],[12,14],[14,25],[25,np.inf]]

AEP = computeAEP(Umean,Ustd,reliability,pwrs,wsps)

print('AEP = {:.4f}kWh'.format(AEP))