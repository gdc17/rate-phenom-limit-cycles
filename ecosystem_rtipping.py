#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 21:02:10 2021

@author: gdc17
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as lin
import imageio

rate = 0.11

def ecosystem_prm(rate):
    
    r1 = 0.95
    r2 = 1.16
    def r(t):
        return r1 + (r2-r1)*(1 + np.tanh(rate*t))/2
    
    m1 = 0.215
    m2 = 0.195
    def m(t):
        return  m1 + (m2-m1)*(1 + np.tanh(rate*t))/2
    
    #------------------- Set up ---------------------------------------------------
    dt = 0.02
    T = 900
    N = int(T/dt)
    K = 100
    Kd = 10000
    time = np.linspace(-T/2, T/2, num=N)
    dt_samp = 40
    i_minus = 0
    i_plus = N
    N_samp = int(N/dt_samp)
    P_vid = np.zeros((N_samp, Kd))
    H_vid = np.zeros((N_samp, Kd))
    
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
    
    def f_P(P,H,r):
        return r*P - C*P**2 - H*g(P)
    
    def f_H(P,H,m):
        return (E*np.exp(-b*P)*g(P)-m)*H   
    
    #------------------ Solve deterministic systen forward in time  ---------------
    P0 = np.random.uniform(15,22, size = Kd); P = P0
    H0 = np.random.uniform(12,18, size = Kd); H = H0
    for i in range(N-1):
         t = time[i]
         P_o = P;
         H_o = H;
         P = P_o + dt*f_P(P_o, H_o, r(t))
         H = H_o + dt*f_H(P_o, H_o, m(t))
         if np.mod(i-i_minus, dt_samp)==0:
            j = int(i/dt_samp)
            P_vid[j,:] = P
            H_vid[j,:] = H
         
    plt.scatter(P, H, s=.3, color='black')
    final = np.stack([P, H]).T
    initial = np.stack([P0, H0]).T
    
    #------------ Circle map ------------------------------------------------------
    
    init_ind = 400
    init = np.stack([P_vid[400,:], H_vid[400]]).T
    c_init = np.array([15,17])
    init_c = init - c_init
    init_thet = np.arctan2(init_c[:,1], init_c[:,0])
    
    for k in range(Kd):
        if final[k,1] < 5:
            final[k] = np.nan
    
    c_final = np.array([15,26])
    final_c = final - c_final
    
    def circ_diff(x,y):
        d1 = x - y
        d2 = -(y - x + 2*np.pi)
        d3 = x - y + 2*np.pi
        if np.abs(d1) <= np.abs(d2) and np.abs(d1) <= np.abs(d3):
            return d1
        elif np.abs(d2) <= np.abs(d1) and np.abs(d2) <= np.abs(d3):
            return d2
        else: 
            return d3
    
    final_thet = np.arctan2(final_c[:,1], final_c[:,0])
    sort_inds = np.argsort(init_thet)
    theta0 = init_thet[sort_inds]
    theta1 = final_thet[sort_inds]
    
    dF = np.zeros(Kd-1)
    dtheta = np.diff(theta0)
    for i in range(Kd-1):
        dF[i] = circ_diff(theta1[i+1], theta1[i])
            
    F_dot = dF/dtheta
    deg = np.round(np.sum(dF)/(2*np.pi))
    length = np.nansum(np.abs(F_dot)*dtheta)/(2*np.pi)
    tip = np.any(np.isnan(theta1))

    #plt.figure()
    #plt.scatter(init_thet, final_thet, s=0.01)
    #plt.xlim([-np.pi,np.pi])
    
    return [theta0, theta1, deg, length, tip, F_dot]


#----------- Phase sensitivity plot -------------------------------------------
K = 200
L = np.zeros(K)
T = np.zeros(K)
r = np.linspace(0.03, 0.12, K)
thet_set = np.zeros((K,3,Kd-1))

for i in range(K):
    print(i)
    theta0, theta1, deg, L[i], T[i], F_dot = ecosystem_prm(r[i])
    thet_set[i,0,:] = r[i]
    thet_set[i,1,:] = theta0[1:]
    #F_dot[F_dot > 40] = 40;
    thet_set[i,2,:] = np.abs(F_dot)
    if L[i] >1 and L[i-1] < 1:
        r1 = i
    if T[i] == True and T[i-1] == False:   
        r2 = i
    
r_scat = np.ndarray.flatten(thet_set[:,0,:])
thet0_scat = np.ndarray.flatten(thet_set[:,1,:])
thet1_scat = np.ndarray.flatten(thet_set[:,2,:])
thet1_scat[thet1_scat > 1000] = 1000

plt.figure(figsize=(10,4))
plt.scatter(r_scat, thet0_scat, c=np.log(thet1_scat), marker='o', s=1)
plt.colorbar(label='$\log|\\frac{dF_r}{d\\theta}$| (Phase Sensitivity)')
plt.xlabel('$r$'); plt.ylabel('$\\theta$')
plt.text(0.1, -2, 'Dangerous Tipping')
plt.axvline(x=0.058, linestyle='--', color='black')
plt.text(0.059, 1.5, 'Rate induced Phase Sensitivity', rotation ='vertical')

#--------------- Video --------------------------------------------------------

x = P_vid; y = H_vid

video = 1
if video == 1:
    # Set spatial resolution
    x_res = 600
    y_res = 400
    xedges = np.linspace(0,60, num=x_res+1)
    yedges = np.linspace(0,40, num=y_res+1)
    # Convert to histogram like data
    H  = np.zeros((len(x), y_res,x_res, 3))
    for i in range(len(x)):
        print(i)
        H[i,:,:,1] = np.histogram2d(y[i,:], x[i,:], bins=(yedges, xedges))[0]
        
    H[H>0] = 1 ; H = H*255; 
    H=H.astype(dtype='uint8')
    imageio.mimwrite('ecosystem_rtipping.mp4', H , fps = 25)

#---------------- Plotting the pullback attractor -----------------------------
T1 = 300; T2 = 550; T3 = 600; T4 = 650; T5 = 700; T6 = 900
plt.scatter(P_vid[T1,:], H_vid[T1,:], s=0.05, label='$t=$')
plt.scatter(P_vid[T2,:], H_vid[T2,:], s=0.05, label='$t=$')
plt.scatter(P_vid[T3,:], H_vid[T3,:], s=0.05, label='$t=$')
plt.scatter(P_vid[T4,:], H_vid[T4,:], s=0.05, label='$t=$')
plt.scatter(P_vid[T5,:], H_vid[T5,:], s=0.05, label='$t=$')
plt.scatter(P_vid[T6,:], H_vid[T6,:], s=0.05, label='$t=$')

plt.xlabel('$P$  (Plant Biomass)'); plt.ylabel('$H$ (Herbivore Biomass)')
plt.ylim([0,35]); plt.xlim([0,60])

