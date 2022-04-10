#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 20:12:16 2022

@author: gdc17
"""
import limit_cycle_functions as lc 
import imageio

dt = 0.001      # time-step length
N = 20000       # nunmber of time-steps
K = 50000       # number of trajectories in ensemble 
dt_samp = 50  # sampling frequency of output 

bet = 2/3       # controls position of repelling limit cycle
b = 2           # shear parameter
r = 1.5         # rate of parameter shift

name = 'pba_r3'  # file name 

Z = lc.pba(dt, N, K, dt_samp, r, bet, b)
x = Z[0,:,:]; y = Z[1,:,:]

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
imageio.mimwrite('%s..mp4' %name, H , fps = 25)