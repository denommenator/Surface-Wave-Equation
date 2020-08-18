#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 21:26:11 2020

@author: robertdenomme
"""

import scipy.ndimage as ndimage
import numpy as np
import matplotlib.pyplot as plt


N=3


K = np.array(
    [[0,1,0],
    [0,0,0],
    [0,1,0]])

Id = np.identity(N)





Ones = ndimage.convolve(Id,K, mode='constant')



L = Ones - 4*Id


L2 = -4*np.identity(N*N)

for i in range(N*N):
    for j in range(N*N):
        if i-j in {1,-1,N,-N}:
            L2[i][j] = 1

Lambda = np.linalg.eig(L)[0].real

Lambda2 = np.linalg.eig(L2)[0].real

fig, ax = plt.subplots()

Lambda2_prime=[]
for l1 in Lambda:
    for l2 in Lambda:
        Lambda2_prime.append(l1+l2)
        
plt.scatter(Lambda2, 0*Lambda2)

plt.scatter(Lambda2_prime, 0*Lambda2+1)

