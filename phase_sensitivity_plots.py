#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 17:40:56 2021

@author: gdc17
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio
import limit_cycle_functions as lc 


dt = 0.001      # time-step length
N = 20000       # nunmber of time-steps
K = 3000         # number of trajectories in ensemble 
dt_samp = N     # we only need the first and last fibres 
bet = 2/3       # controls position of repelling limit cycle
b = 2           # shear parameter
    
#------- Length plots ---------------------------------------------------------
K_R = 300                        # number of r values to use 
L = np.zeros(K_R)                # stores average phase sensitivity 
T = np.zeros(K_R)                # stores Boolean tip or not 
r = np.linspace(0.5, 3, K_R)     # array of r values 
thet_set = np.zeros((K,3,K-1))


# loop through r values and store phase response map data 
for i in range(K_R):
    print(i)
    Z = Z = lc.pba(dt, N, K, dt_samp, r[i], bet, b)
    theta0, theta1, deg, L[i], T[i], F_dot = lc.prm(Z, extra = 1)
    thet_set[i,0,:] = r[i]
    thet_set[i,1,:] = theta0[1:]
    #F_dot[F_dot > 40] = 40;
    thet_set[i,2,:] = np.abs(F_dot)
    if L[i] >1 and L[i-1] < 1:
        r1 = i
    if T[i] == True and T[i-1] == False:
        r2 = i
    
# convert data into form for scatter plot 
r_scat = np.ndarray.flatten(thet_set[:,0,:])
thet0_scat = np.ndarray.flatten(thet_set[:,1,:])
thet1_scat = np.ndarray.flatten(thet_set[:,2,:])
thet1_scat[thet1_scat > 500] = 500 # caps out the infinite phase sensitivty 
thet0_scat = thet0_scat + np.pi

# plot r against the average phase sensitivity (length of phase response map)
plt.figure(figsize=(10,4))
plt.scatter(r_scat, thet0_scat, c=np.log(thet1_scat), marker='o', s=.1)
plt.colorbar(label='$\log|\\frac{dF_r}{d\\theta}$| (Log Phase Sensitivity)')
plt.xlabel('$r$'); plt.ylabel('$\\theta$')
# labels are currently harcoded in
plt.text(2.2, 0.1, 'Dangerous Tipping')
plt.text(2.5, 2.7, 'Safe Tipping')
plt.arrow(2.64,2.66, -0.08, -0.45, head_width=.04)
plt.axvline(x=0.88, linestyle='--', color='black')
plt.text(0.81, 1.7+np.pi, 'Phase Sensitivity Bifurcation', rotation ='vertical')
plt.show()

# produce plot showing phase sensitivity of each phase
plt.figure(figsize=(10,3.5))
plt.axvline(x=0.88, linestyle='--', color='blue')
plt.text(0.85, 7.2, 'Phase Sensitivity Bifurcation', rotation ='vertical', color='blue')
plt.axvline(x=r[r2], linestyle='--', color='red')
plt.text(1.48, 4.5, 'Tipping', rotation ='vertical', color='red')
#plt.fill_between(r[r1-1:r2+1], 0,10, alpha=.5, color='orange', label='Stretching and Folding')
#plt.fill_between(r[r2:K], 0,10, alpha=.5, color='red', label='Tipping')
#plt.fill_between(r[0:r1], 0,10, alpha=.5, color='yellow', label='Small $r$ Behaviour')
plt.plot(r,L, label='$\\bar{S}_r$', color='black', linewidth=2)
plt.ylim([0,8]); plt.xlim([0.5,2])
plt.ylabel('$\\bar{S}_r$ (Average Phase Sensitivity)'); plt.xlabel('$r$')
plt.legend()
plt.show()

