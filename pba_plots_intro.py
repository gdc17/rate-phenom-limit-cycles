#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 21:52:41 2022

@author: gdc17
"""

import numpy as np
import matplotlib.pyplot as plt
import limit_cycle_functions as lc

def generate_pba_plots_intro(r, t):
    dt = 0.001     # size of time-step
    N = 50000      # number of time-step
    K = 50000       # number of trajecotires
    dt_samp = 50   # sampling frequency of output 
  #  r = 1.75        # rate of parameter shift
    bet = 2/3      # controls position of repelling limit cycle
    b = 2          # control amount of shear 
    
    # use pba function to obtain approximation of pullback attractor 
    Z = lc.pba(dt, N, K, dt_samp, r, bet, b)
    x = Z[0,:,:]; y = Z[1,:,:]
    
   # t = 560  # time to plot
    time = - 25 + t*dt_samp*dt
    
    # Construct the frozen in time attracting and repelling limit cycle 
    thet = np.linspace(0, 2*np.pi, num=K)
    gamma = np.zeros((2,K)); gamma[0,:] = np.cos(thet); gamma[1,:] = np.sin(thet)
    gamma_r = np.zeros((2,K))
    gamma_r[0,:] = (1/bet)*np.cos(thet); gamma_r[1,:] = (1/bet)*np.sin(thet)
    
    # Plot the t fibre of the pullback attracting along with frozen in time objects
    plt.figure(figsize=(4,4))
    plt.plot(gamma[0,:]+ lc.q(r*time), gamma[1,:], color='blue',
             label='$\Gamma_{\\lambda(rt)}$', linewidth=0.5)
    plt.scatter(x[t,:],y[t,:], s=0.1, color='green')
    plt.plot([1],[1], color='green', label='$A_r(t)$')
    plt.ylim([-1.8, 1.8]), plt.xlim([-1.8+lc.q(r*time), 1.8+lc.q(r*time)])
    plt.text(-1.7+lc.q(r*time), -1.7, '$t = $%1.1f' %time, fontsize=14 )
    plt.text(-1.7+lc.q(r*time), -1.3, '$r = $%1.2f' %r, fontsize=14 )
    plt.legend(loc='upper right'); plt.xlabel('$x_1$'); plt.ylabel('$x_2$')
    plt.show()

def plot():
    r_vals = [0.8, 1.25, 1.75]
    t_vals = [480,520,560]
    for r in r_vals:
        for t in t_vals:
            generate_pba_plots_intro(r, t)
   
if __name__ == '__main__':  
    plot()


