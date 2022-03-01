#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:58:54 2021

@author: gdc17
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as lin

r = 1
m = 0.21

#------------------- Set up ---------------------------------------------------
dt = 0.01
T = 400
N = int(T/dt)
K = 100
Kd = 800000

#---------------- Parameters --------------------------------------------------
C = 0.02
a = 10
b =  0.005
b_c =  0.01
E = 0.4
c_max = 1

#--------------- Define functions ---------------------------------------------
def g(P):
    return c_max*(P**2/(P**2+a**2))*np.exp(-b_c*P)

def Dg(P):
    return -g(P)*(2*P/(P**2+a**2)-2/P+b_c)

def f_P(P,H,r):
    return r*P - C*P**2 - H*g(P)

def f_H(P,H,m):
    return (E*np.exp(-b*P)*g(P)-m)*H
 
def Df(P,H,r,m):
    a = r - 2*C*P - H*Dg(P)
    b = -g(P)
    c = H*E*np.exp(-b*P)*(Dg(P)-b*g(P))
    d = E*g(P)*np.exp(-b*P) - m
    return np.matrix([[a,b],[c,d]])

#------------------ Solve deterministic systen forward in time  ---------------
P0 = np.random.uniform(20,21, size = Kd); P = P0
H0 = np.random.uniform(12,14, size = Kd); H = H0
for i in range(N-1):
     t = i*dt
     P_o = P;
     H_o = H;
     P = P_o + dt*f_P(P_o, H_o, r)
     H = H_o + dt*f_H(P_o, H_o, m)
     
plt.scatter(P, H, s=.3, color='black')
final = np.stack([P, H]).T
initial = np.stack([P0, H0]).T

l1 = np.random.randint(0,Kd)
l2 = np.random.randint(0,Kd)
l3 = np.random.randint(0,Kd)
l4 = np.random.randint(0,Kd)
l5 = np.random.randint(0,Kd)
l6 = np.random.randint(0,Kd)
eps = .3
Ws_ind1 = []
Ws_ind2 = []
Ws_ind3 = []
Ws_ind4 = []
Ws_ind5 = []
Ws_ind6 = []
for k in range(Kd):
    if np.sum(np.abs(final[l1]-final[k])) < eps:
        Ws_ind1.append(k)
    if np.sum(np.abs(final[l2]-final[k])) < eps:
        Ws_ind2.append(k)
    if np.sum(np.abs(final[l3]-final[k])) < eps:
        Ws_ind3.append(k)
    if np.sum(np.abs(final[l4]-final[k])) < eps:
        Ws_ind4.append(k)
    if np.sum(np.abs(final[l5]-final[k])) < eps:
        Ws_ind5.append(k)
    if np.sum(np.abs(final[l6]-final[k])) < eps:
        Ws_ind6.append(k)
        
Ws1 = initial[Ws_ind1]
Ws2 = initial[Ws_ind2]
Ws3 = initial[Ws_ind3]
Ws4 = initial[Ws_ind4]
Ws5 = initial[Ws_ind5]
Ws6 = initial[Ws_ind6]

plt.figure(figsize=(8,8))
plt.scatter(P, H, s=.3, color='black')
plt.scatter(Ws1[:,0], Ws1[:,1], s=.3, color='darkorange')   
plt.scatter(Ws2[:,0], Ws2[:,1], s=.3, color='blue')      
plt.scatter(Ws3[:,0], Ws3[:,1], s=.3, color='green')  
plt.scatter(Ws4[:,0], Ws4[:,1], s=.3, color='magenta')  
plt.scatter(Ws5[:,0], Ws5[:,1], s=.3, color='greenyellow') 
plt.scatter(Ws6[:,0], Ws6[:,1], s=.3, color='aqua') 

plt.plot(Pb, Hb, linewidth = 1, color = 'red')
plt.text(20.4,12.25, '$W^s(x_s)$', fontsize=15, color='red')
plt.ylim([12,14]); plt.xlim([20,21])
plt.xlabel('$P$  (Plant Biomass)'); plt.ylabel('$H$ (Herbivore Biomass)')

