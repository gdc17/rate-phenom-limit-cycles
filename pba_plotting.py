#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 20:10:49 2022

@author: gdc17
"""

 
        
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


