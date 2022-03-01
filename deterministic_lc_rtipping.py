#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:13:03 2021

@author: gdc17
"""

import numpy as np
import matplotlib.pyplot as plt
import imageio

# Time parameters
dt = 0.001
N = 50000
K = 8000
T = N*dt
time = np.linspace(-T/2,T/2,N)
i_minus = int(N/2) - 25000
i_plus =  int(N/2) + 25000
dt_samp = 50
N_samp = int((i_plus - i_minus)/(dt_samp))

# Parameters
r = 1.75; C = 2
bet = 2/3; a = 1
b = 2; 

# Initialise arrays
Z = np.zeros((2,K))
r0 = np.random.uniform(0.9, 1.1, size=K)
thet = np.random.uniform(0, 2*np.pi, size=K)
Z[0,:] = r0*np.cos(thet); Z[1,:] = r0*np.sin(thet)
Z_vid = np.zeros((2, N_samp, K))

 
def q(t):                                 
   return C*(np.tanh(r*t)+ 1)

def f(Z):
    x = Z[0]; y = Z[1]
    r2 = x**2 + y**2
    r = np.sqrt(r2)
    thet = np.arctan2(y,x)
    f1 = x*(1 - r)*(1 - bet*r) - y*(1 + b*r2)
    f2 = y*(1 - r)*(1 - bet*r) + x*(1 + b*r2)
    return np.array([f1,f2])

#----------------- Solve system ------------------------------------------        
for i in range(N-1):
    if np.random.normal() > 2:
        print(i)
    qi = q(time[i])*np.ones(K)
    qix =  np.stack((qi, np.zeros(K)), 1).T
    Z_o = Z;
    Z = Z_o + dt*f(Z_o - qix)  
    # save data for video 
    if i >= i_minus and i < i_plus and np.mod(i-i_minus, dt_samp)==0:
        j = int((i-i_minus)/dt_samp)
        Z_vid[:,j,:] = Z

#-------------------- Plotting ------------------------------------------------        
x = Z_vid[0,:,:]; y = Z_vid[1,:,:]

# Plot pullback attractor and repellor
t = 560
time = - 25 + t*dt_samp*dt
thet = np.linspace(0, 2*np.pi, num=K)
gamma = np.zeros((2,K))
gamma_r = np.zeros((2,K))
gamma[0,:] = np.cos(thet); gamma[1,:] = np.sin(thet)
gamma_r[0,:] = (1/bet)*np.cos(thet); gamma_r[1,:] = (1/bet)*np.sin(thet)


plt.figure(figsize=(4,4))
plt.plot(gamma[0,:]+ q(time), gamma[1,:], color='blue', label='$\Gamma_{\\lambda(rt)}$', linewidth=0.5)
plt.plot(gamma_r[0,:]+q(time), gamma_r[1,:], color='red', label='$\\Upsilon_{\\lambda(rt)}$', linewidth=0.5)
plt.scatter([q(time)], [0], color='black', label='$R_{\lambda(rt)}$')
plt.scatter(x[t,:],y[t,:], s=0.1, color='green')
plt.plot([1],[1], color='green', label='$A_r(t)$')
plt.ylim([-1.8, 1.8]), plt.xlim([-1.8+q(time), 1.8+q(time)])
plt.text(-1.7+q(time), -1.7, '$t = $%1.1f' %time, fontsize=14 )
plt.text(-1.7+q(time), -1.3, '$r = $%1.2f' %r, fontsize=14 )
plt.legend(loc='upper right')

#----------------- Video ------------------------------------------------------

video = 1
if video == 1:
    # Set spatial resolution
    x_res = 600
    y_res = 425
    xedges = np.linspace(-1.7,5.7, num=x_res+1)
    yedges = np.linspace(-2,2, num=y_res+1)
    # Convert to histogram like data
    H  = np.zeros((len(x), y_res,x_res, 3))
    for i in range(len(x)):
        print(i)
        H[i,:,:,1] = np.histogram2d(y[i,:], x[i,:], bins=(yedges, xedges))[0]
        
    H[H>0] = 1 ; H = H*255
    H=H.astype(dtype='uint8')
    imageio.mimwrite('pba_r3.mp4', H , fps = 25)