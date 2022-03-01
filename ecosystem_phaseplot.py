#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 18:13:33 2021

@author: gdc17
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 20:58:54 2021

@author: gdc17
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as lin
import scipy.optimize as op

r = 1
m = 0.21

#------------------- Set up ---------------------------------------------------
dt = 0.01
T = 500
N = int(T/dt)
K1 = 100
K2 = 500

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

#-------------- Approximate equilibra and stable/unstable manifold------------
def F(x):
    x0 = x[0]; x1 = x[1]
    F1 = f_P(x0, x1, r)
    F2 = f_H(x0, x1, m)
    return np.array([F1,F2])

xs = op.fsolve(F, [40,15])
xex = op.fsolve(F, [50,0])

N = int(T/dt)
P = np.zeros((N,K1)); Pb = np.zeros((N,K1))
H = np.zeros((N,K1)); Hb = np.zeros((N,K1))

eps = 0.01
P[0,:] = np.random.uniform(xs[0]-eps,  xs[0]+eps, size = K1); 
Pb[0,:] = np.random.uniform(xs[0]-eps, xs[0]+eps, size = K1); 
H[0,:] = np.random.uniform(xs[1]-eps,  xs[1]+eps, size = K1); 
Hb[0,:] = np.random.uniform(xs[1]-eps, xs[1]+eps, size = K1); 
for i in range(N-1):
     t = i*dt
     P[i+1,:] = P[i,:] + dt*f_P(P[i,:], H[i,:], r) 
     Pb[i+1,:] = Pb[i,:] - dt*f_P(Pb[i,:], Hb[i,:], r) 
     H[i+1,:] = H[i,:] + dt*f_H(P[i,:], H[i,:], m)
     Hb[i+1,:] = Hb[i,:] - dt*f_H(Pb[i,:], Hb[i,:], m) 
     
#---------- Plot other trajectories -------------------------------------------
P1 = np.zeros((N,K2)); P1b = np.zeros((N,K2))
H1 = np.zeros((N,K2)); H1b = np.zeros((N,K2))     

P1[0,:] = np.random.uniform(0,  52, size = K2); 
H1[0,:] = np.random.uniform(0,  28, size = K2); 
for i in range(N-1):
     t = i*dt
     P1[i+1,:] = P1[i,:] + dt*f_P(P1[i,:], H1[i,:], r) 
     H1[i+1,:] = H1[i,:] + dt*f_H(P1[i,:], H1[i,:], m)

plt.figure(figsize=(10,6))
plt.plot(P1,H1, linewidth = .1)
plt.plot(P,H, linewidth = 1, color = 'blue')
plt.plot(Pb, Hb, linewidth = 1, color = 'red')
plt.scatter(xs[0], xs[1], color = 'black')
plt.scatter(xex[0], xex[1], color = 'black')
plt.ylim([0,28]); plt.xlim([0,52])
plt.text(15,10.5, '$W^s(x_s)$', fontsize=15, color='red')
plt.text(30.7,20, '$W^u(x_s)$', fontsize=15, color='blue')
plt.text(39,14.5, '$x_s$', fontsize=15, color='black')
plt.xlabel('$P$  (Plant Biomass)'); plt.ylabel('$H$ (Herbivore Biomass)')

     




        


