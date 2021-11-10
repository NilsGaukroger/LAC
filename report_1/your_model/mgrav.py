# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 10:10:04 2021

@author: nilsg
"""

import pandas as pd
import numpy as np

filepath = "./results_redesign/info/inertia.dat"
g = 9.81

#%% Import data

rows_to_skip = [i for i in range(8)]
rows_to_skip.append(8)
rows_to_skip.extend([i for i in range(18,39)])

df = pd.read_csv(filepath,names=["Body name", "mass", "cog_global_x", "cog_global_y", "cog_global_z", "cog_local_x", "cog_local_y", "cog_local_z"],skiprows=rows_to_skip,delim_whitespace=True)

#%% Calculate RNA CoG in global coordinates

df["cog_global"] = np.sqrt(df["cog_global_x"]**2 + df["cog_global_y"]**2 + df["cog_global_z"]**2)

CoG_computed = [((df["mass"] * df["cog_global_x"]) / df["mass"].sum()).sum()]
CoG_computed.append(((df["mass"] * df["cog_global_y"]) / df["mass"].sum()).sum())
CoG_computed.append(((df["mass"] * df["cog_global_z"]) / df["mass"].sum()).sum())

#%%

CoG = pd.read_csv(filepath,skiprows=26,nrows=1,header=None,delim_whitespace=True).iloc[:,2:5].values.tolist()[0]

Mgrav = -CoG[1] * df["mass"].sum() * g
print(Mgrav/1e3)