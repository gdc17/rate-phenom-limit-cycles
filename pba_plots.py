#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 21:52:41 2022

@author: gdc17
"""

import numpy as np
import matplotlib.pyplot as plt
import limit_cycle_functions as lc

dt = 0.001     # size of time-step
N = 50000      # number of time-step
K = 8000       # number of trajecotires
dt_samp = 50   # sampling frequency of output 
r = 1.4        # rate of parameter shift
bet = 2/3      # controls position of repelling limit cycle
b = 2          # control amount of shear 

# use pba function to obtain approximation of pullback attractor 
Z = lc.pba(dt, N, K, dt_samp, r, bet, b)
x = Z[0,:,:]; y = Z[1,:,:]

t = 570  # time to plot
time = - 25 + t*dt_samp*dt

# Construct the frozen in time attracting and repelling limit cycle 
thet = np.linspace(0, 2*np.pi, num=K)
gamma = np.zeros((2,K)); gamma[0,:] = np.cos(thet); gamma[1,:] = np.sin(thet)
gamma_r = np.zeros((2,K))
gamma_r[0,:] = (1/bet)*np.cos(thet); gamma_r[1,:] = (1/bet)*np.sin(thet)

# Plot the t fibre of the pullback attracting along with frozen in time objects
plt.figure(figsize=(4,4))
plt.plot(gamma[0,:]+ q(time), gamma[1,:], color='blue',
         label='$\Gamma_{\\lambda(rt)}$', linewidth=0.5)
plt.plot(gamma_r[0,:]+q(time), gamma_r[1,:], color='red',
         label='$\\Upsilon_{\\lambda(rt)}$', linewidth=0.5)
plt.scatter([q(time)], [0], color='black', label='$R_{\lambda(rt)}$')
plt.scatter(x[t,:],y[t,:], s=0.1, color='green')
plt.plot([1],[1], color='green', label='$A_r(t)$')
plt.ylim([-1.8, 1.8]), plt.xlim([-1.8+q(time), 1.8+q(time)])
plt.text(-1.7+q(time), -1.7, '$t = $%1.1f' %time, fontsize=14 )
plt.text(-1.7+q(time), -1.3, '$r = $%1.2f' %r, fontsize=14 )
plt.legend(loc='upper right')
plt.show()



