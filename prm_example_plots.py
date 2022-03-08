#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 19:36:27 2022

@author: gdc17
"""

import numpy as np
import matplotlib.pyplot as plt
import limit_cycle_functions as lc

dt = 0.001      # time-step length
N = 20000       # nunmber of time-steps
K = 50000       # number of trajectories in ensemble 
dt_samp = 1000  # sampling frequency of output 

bet = 2/3       # controls position of repelling limit cycle
b = 2           # shear parameter
r = 1.5         # rate of parameter shift


Z = lc.pba(dt, N, K, dt_samp, r, bet, b)
theta0, theta1 = lc.prm(Z, extra = 0 )
plt.figure(figsize=(4,4))
plt.scatter(theta0, theta1, s=0.001)
plt.xlabel('$\\theta$'); plt.ylabel('$F_r(\\theta)$')
plt.show()
