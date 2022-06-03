#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 16:38:24 2021

@author: gdc17
"""
beta = 1.5
b = -1.6
M = 400000
r = np.linspace(0,beta, num=M)
def g(r):
    return b*(np.log(r)/beta - np.log(beta-r)*(beta+1)/beta)

thet0 = 0
iso0 = np.zeros((2,M))
iso0[0,:] = r*np.cos(thet0 + g(r))
iso0[1,:] = r*np.sin(thet0 + g(r))

thet1 = 2*np.pi/6
iso1 = np.zeros((2,M))
iso1[0,:] = r*np.cos(thet1 + g(r))
iso1[1,:] = r*np.sin(thet1 + g(r))

thet2 = 4*np.pi/6
iso2 = np.zeros((2,M))
iso2[0,:] = r*np.cos(thet2 + g(r))
iso2[1,:] = r*np.sin(thet2 + g(r))

thet3 = 6*np.pi/6
iso3 = np.zeros((2,M))
iso3[0,:] = r*np.cos(thet3 + g(r))
iso3[1,:] = r*np.sin(thet3 + g(r))

thet4 = 8*np.pi/6
iso4 = np.zeros((2,M))
iso4[0,:] = r*np.cos(thet4 + g(r))
iso4[1,:] = r*np.sin(thet4 + g(r))

thet5 = 10*np.pi/6
iso5 = np.zeros((2,M))
iso5[0,:] = r*np.cos(thet5 + g(r))
iso5[1,:] = r*np.sin(thet5 + g(r))


C = 200
plt.figure(figsize=(8,8))
plt.plot(iso0[0,:-C], iso0[1,:-C], linewidth=1, color='darkorange')   
plt.plot(iso1[0,:-C], iso1[1,:-C], linewidth=1, color='blue')
plt.plot(iso2[0,:-C], iso2[1,:-C], linewidth=1, color='green')
plt.plot(iso3[0,:-C], iso3[1,:-C], linewidth=1, color='magenta')
plt.plot(iso4[0,:-C], iso4[1,:-C], linewidth=1, color='greenyellow')
plt.plot(iso5[0,:-C], iso5[1,:-C], linewidth=1, color = 'aqua')


K = 10000
thet = np.linspace(0, 2*np.pi, num=K)
gamma = np.zeros((2,K))
gamma_r = np.zeros((2,K))
gamma[0,:] = np.cos(thet); gamma[1,:] = np.sin(thet)
gamma_r[0,:] = beta*np.cos(thet); gamma_r[1,:] = beta*np.sin(thet)

plt.plot(gamma[0,:], gamma[1,:], linewidth=2,color='black', label='$\Gamma$')
plt.plot(gamma_r[0,:], gamma_r[1,:], linewidth=3, color='red', label='$W^s(\Gamma)$')
plt.ylabel('$x_2$', fontsize=14); plt.xlabel('$x_1$', fontsize=14)
plt.legend()
#plt.ylim([-1.52, -1.4]); plt.xlim([-0.1, 0.1])


