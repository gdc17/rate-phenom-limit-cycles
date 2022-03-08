#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 19:58:21 2022

@author: gdc17
"""


import numpy as np
import matplotlib.pyplot as plt


def q(t):
    ''' Parameter shift function '''                                 
    return 2*(np.tanh(r*t)+ 1)


def pba(dt = 0.01, N = 50000, K = 8000, dt_samp = 50,
        r = 1.75, bet = 2/3, b = 2):
    '''Returns a numerical approximation of an ensemble of K randomly chosen
       trajectories on the pullback attractor '''

    # Time parameters
    T = N*dt
    time = np.linspace(-T/2, T/2, N)
    if N % dt_samp != 0:
         raise ValueError("Sampling frequency should divide number of time-steps")
    N_samp = int(N/dt_samp + 1)
    
    # Initialise arrays
    Z = np.zeros((2,K))
    r0 = np.random.uniform(1, 1, size = K)
    thet = np.random.uniform(0, 2*np.pi, size=K)
    Z[0,:] = r0*np.cos(thet); Z[1,:] = r0*np.sin(thet)
    Z_vid = np.zeros((2, N_samp, K))
    
    # define parameter shift function
    def q(t):                                 
       return 2*(np.tanh(r*t)+ 1)
    
    # define vector field (to be shifted)
    def f(Z):
        x = Z[0]; y = Z[1]
        r2 = x**2 + y**2
        r = np.sqrt(r2)
        f1 = x*(1 - r)*(1 - bet*r) - y*(1 + b*r2)
        f2 = y*(1 - r)*(1 - bet*r) + x*(1 + b*r2)
        return np.array([f1,f2])
    
    # Solve system (using Euler mathod) 
    j = 0 
    for i in range(N-1):
        qi = q(time[i])*np.ones(K)
        qi_stacked =  np.stack((qi, np.zeros(K)), 1).T
        Z_o = Z;
        Z = Z_o + dt*f(Z_o - qi_stacked)  
        # save data for video 
        if np.mod(i, dt_samp) == 0 or i == N-2: 
            Z_vid[:,j,:] = Z
            j = j + 1
        
    return Z_vid

#------------------------------------------------------------------------------

def prm(pba, extra = 1):  
    '''Returns the phase response map and takes the output of the pba function
       as an input. Returns the data of the phase response along with some
       extra properties'''
    
    # convert from cartesian to angular variables 
    x = pba[0,:,:]; y = pba[1,:,:]    
    K = np.shape(pba)[2]
    theta0 = np.arctan2(y[0,:], x[0,:])
    theta1 = np.arctan2(y[-1,:], x[-1,:]-2*2)
    sort_inds = np.argsort(theta0)
    theta0 = theta0[sort_inds]
    theta1 = theta1[sort_inds]
    
    def circ_diff(x,y):
        ''' Defines flat metric on circle parametersied from 0 to 2pi' to 
            calculate derivatives'''
        d1 = x - y
        d2 = -(y - x + 2*np.pi)
        d3 = x - y + 2*np.pi
        if np.abs(d1) <= np.abs(d2) and np.abs(d1) <= np.abs(d3):
            return d1
        elif np.abs(d2) <= np.abs(d1) and np.abs(d2) <= np.abs(d3):
            return d2
        else: 
            return d3
     
    if extra == 1:
        dF = np.zeros(K-1)
        lift = np.zeros(K)
        dtheta = np.diff(theta0)
        for i in range(K-1):
            dF[i] = circ_diff(theta1[i+1], theta1[i])
            lift[i+1] = lift[i] + dF[i]
        F_dot = dF/dtheta
        deg = np.round(np.sum(dF)/(2*np.pi))
        length = np.nansum(np.abs(F_dot)*dtheta)/(2*np.pi)
        tip = np.any(np.isnan(theta1))
        return [theta0, theta1, deg, length, tip, F_dot]
    else:
        return [theta0,theta1]
    