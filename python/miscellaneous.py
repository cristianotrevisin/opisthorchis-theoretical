#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:11:08 2023

@author: cristiano
"""

import numpy as np
import matplotlib.pyplot as plt

def zipf(rank, expn, minP):
    H = np.sum(1. / (rank ** expn))
    P = 1. / (rank ** expn) / H
    P = P / np.min(P) * minP
    return P

def drawborders(MSC):
    for ix in range(1, MSC.shape[0]):
        for iy in range(1, MSC.shape[1]):
            if MSC[ix, iy] != MSC[ix, iy-1]:
                plt.plot([ix-0.5, ix+0.5], [iy-0.5, iy-0.5], color='black', linewidth=0.5)
                
    for iy in range(1, MSC.shape[1]):
        for ix in range(1, MSC.shape[0]):
            if MSC[ix, iy] != MSC[ix-1, iy]:
                plt.plot([ix-0.5, ix-0.5], [iy-0.5, iy+0.5], color='black', linewidth=0.5)
                

def assignSC(MSC, X):
    M = np.zeros_like(MSC)
    for iz in range(len(X)):
        M[MSC == iz] = X[iz]
    return M
