#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 13:58:18 2025

@author: hmonteiro
"""

import numpy as np
import matplotlib.pyplot as plt


h2ofile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/h20_crystalline.nk'
h2ofile08 = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/h20_Warren08.nk'



h2Odata = np.loadtxt(h2ofile,max_rows=2336,skiprows=2)
h2O08data = np.loadtxt(h2ofile08,max_rows=486,skiprows=2)

plt.figure()
plt.xscale('log')
# plt.yscale('log')

plt.plot(h2Odata[:,0], h2Odata[:,1],'C0',label='n')
plt.plot(h2Odata[:,0], h2Odata[:,2],'C1',label='k')

plt.plot(h2O08data[:,0], h2O08data[:,1],'C0--',label='n')
plt.plot(h2O08data[:,0], h2O08data[:,2],'C1--',label='k')

newlambda = np.concatenate((h2Odata[h2Odata[:,0]< 0.0443,0], h2O08data[:,0]))
new_n = np.concatenate((h2Odata[h2Odata[:,0]< 0.0443,1], h2O08data[:,1]))
new_k = np.concatenate((h2Odata[h2Odata[:,0]< 0.0443,2], h2O08data[:,2]))

plt.plot(newlambda, new_n,'C2',label='n')
plt.plot(newlambda, new_k,'C2',label='k')


header = "nk\nh20_ice.nk 273 0.917 0.588 6.005"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/H2O_ice.nk'
np.savetxt(output_file, np.vstack((newlambda, new_n, new_k)).T, header=header, comments='', fmt='%.6E')



















