#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 17:40:56 2021

@author: gdc17
"""
import numpy as np
import matplotlib.pyplot as plt
import imageio

def tipping_circle_map(r,b):
    ### Time parameters
    dt = 0.001
    N = 20000
    K = 30000
    T = N*dt
    time = np.linspace(-T/2,T/2,N)
    dt_samp = 1000
    N_samp = int(N/dt_samp)
    
    # Parameters
    r = r; C = 2
    bet = 2/3; b = b
    
    # Initialise arrays
    Z = np.zeros((2,K))
    r0 = np.random.uniform(1,1, size=K)
    thet = np.random.uniform(0*np.pi,2*np.pi, size=K)
    Z[0,:] = r0*np.cos(thet)
    Z[1,:] = r0*np.sin(thet)
    Z_vid = np.zeros((2, N_samp, K))
        
    def q(t):                                 
       return C*(np.tanh(r*t) + 1)
        
    def f(Z):
        x = Z[0]; y = Z[1]
        r2 = x**2 + y**2
        r = np.sqrt(r2)
        thet = np.arctan2(y,x)
        f1 = x*(1-r)*(1-bet*r) - y*(1+b*r2)
        f2 = y*(1-r)*(1-bet*r) + x*(1+b*r2)
        return np.array([f1,f2])
    
    #----------------- Solve system -------------------------------------------        
    for i in range(N-1):
        qi = q(time[i])*np.ones(K)
        qix =  np.stack((qi, np.zeros(K)), 1).T
        Z_o = Z;
        Z = Z_o + dt*f(Z_o - qix)  
        # save data for video 
        if np.mod(i-i_minus, dt_samp)==0:
            j = int(i/dt_samp)
            Z_vid[:,j,:] = Z
    
    #-------------------- Plotting --------------------------------------------        
    x = Z_vid[0,:,:]; y = Z_vid[1,:,:]    

    #plt.figure()
    theta0 = np.arctan2(y[1,:], x[1,:])
    theta1 = np.arctan2(y[-1,:], x[-1,:]-2*C)
    sort_inds = np.argsort(theta0)
    theta0 = theta0[sort_inds]
    theta1 = theta1[sort_inds]
    
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
     
    extra = 1
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
    else:
        deg = np.nan; length = np.nan; tip = np.nan
    

    plot = 1
    if plot == 1:
        plt.scatter(theta0, theta1, s=0.01)
        plt.ylim([-np.pi, np.pi])
        plt.xlim([-np.pi, np.pi])
        
    return [theta0, theta1, deg, length, tip, F_dot]

pre0, pre1 = tipping_circle_map(1.6,-4)[0:2]
post0, post1 = tipping_circle_map(1.7,-4)[0:2]
plt.scatter(pre0, pre1, s=0.01)
plt.scatter(post0, post1, s=0.01)
plt.ylim([-np.pi, np.pi])
plt.xlim([-np.pi, np.pi])

#------- Length plots ---------------------------------------------------------
K = 300
L = np.zeros(K)
T = np.zeros(K)
r = np.linspace(0.5, 3, K)
thet_set = np.zeros((K,3,30000-1))

for i in range(K):
    print(i)
    theta0, theta1, deg, L[i], T[i], F_dot = tipping_circle_map(r[i], 2)
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
thet1_scat[thet1_scat > 500] = 500

thet0_scat = thet0_scat - np.pi
plt.figure(figsize=(10,4))
plt.scatter(r_scat, thet0_scat, c=np.log(thet1_scat), marker='o', s=.1)
plt.colorbar(label='$\log|\\frac{dF_r}{d\\theta}$| (Log Phase Sensitivity)')
plt.xlabel('$r$'); plt.ylabel('$\\theta$')
plt.text(2.2, 0+np.pi, 'Dangerous Tipping')
plt.text(2.5, 2.7+np.pi, 'Safe Tipping')
plt.arrow(2.64,2.66+np.pi, -0.08, -0.45, head_width=.04)
plt.axvline(x=0.88, linestyle='--', color='black')
plt.text(0.81, 1.7+np.pi, 'Phase Sensitivity Bifurcation', rotation ='vertical')

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

#------------- Nice big subplot -----------------------------------------------

fig, axs = plt.subplots(2, gridspec_kw={'height_ratios': [0.7, 2]}, figsize=(10, 6))

plot0 = axs[0].plot(r,L, label='$\\bar{S}_r$ (Average Phase Sensitivity)', color='black', linewidth=2)
axs[0].set_ylim([0,10])
axs[0].set_ylabel('$\\bar{S}_r$')
axs[0].fill_between(r[r1-1:r2+1], 0,10, alpha=.5, color='orange', label='Stretching and Folding')
axs[0].fill_between(r[r2:K], 0,10, alpha=.5, color='red', label='Partial tipping')
axs[0].fill_between(r[0:r1], 0,10, alpha=.5, color='yellow', label='Small $r$ Behaviour')
axs[0].axvline(x=r[r1], linestyle='--', color='black')
axs[0].axvline(x=r[r2], linestyle='--', color='black')
axs[0].legend()

plot1 = axs[1].scatter(r_scat, thet0_scat, c=np.log(thet1_scat), marker='o', s=0.02)
axs[1].set_ylabel('$\\theta$')
axs[1].set_xlabel('$r$')
axs[1].text(2.6, 0.9, 'Dangerous Tipping')
axs[1].axvline(x=1.32, linestyle='--', color='black')
axs[1].text(1.34, 1.4, 'Stretching and Folding', rotation ='vertical')
axs[1].annotate('Safe Tipping (Invisible)', xy=(2.85, -2.4), xytext=(2.5, -1),
                    arrowprops=dict(arrowstyle="->"))
cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.45])
plt.colorbar(plot1, cax=cbar_ax, label='$\log|\\frac{dF_r}{d\\theta}(\\theta)$| (Phase Sensitivity)')

#------------ Phase Response Map ----------------------------------------------
theta0, theta1, deg, length, tip, F_dot = tipping_circle_map(1.75,2)
plt.figure(figsize=(4,4))
plt.scatter(theta0, theta1, s=0.001)
plt.xlabel('$\\theta$'); plt.ylabel('$F_r(\\theta)$')


#-------- Bifurcation heatmaps ------------------------------------------------

#K1 = 25; K2 = 25
#r = np.linspace(0.75, 1.5, K1)
#b = np.linspace(-2,-4,K2)
K1 = 10; K2 = 1
r = np.linspace(1, 3, K1)
b = np.linspace(-2,-4,K2)
A = np.zeros((K1,K2))
B = np.zeros((K1,K2))
C = np.zeros((K1,K2))

for i in range(K1):
    print(i)
    for j in range(K2):
        tcm = tipping_circle_map(r[i],b[j])
        A[i,j] = tcm[2]; B[i,j] = tcm[3]; C[i,j] = tcm[4]

plt.figure()      
plt.imshow(np.flip(A,0), extent=[b[0],b[-1],r[0],r[-1]])
plt.xlabel('$b$'); plt.ylabel('$r$')
plt.figure()
plt.imshow(np.flip(B,0), vmin=0, vmax=6, extent=[b[0],b[-1],r[0],r[-1]])
plt.xlabel('$b$'); plt.ylabel('$r$')
plt.colorbar()
plt.figure()
plt.imshow(np.flip(C,0), extent=[b[0],b[-1],r[0],r[-1]])
plt.xlabel('$b$'); plt.ylabel('$r$')

#------------ Video -----------------------------------------------------------

r = np.linspace(0.57, 0.66, 141)
theta0 = np.zeros((500, len(r)))
theta1 = np.zeros((500, len(r)))


for i in range(len(r)):
    print(i)
    theta0[:,i], theta1[:,i] = tipping_circle_map(r[i], 1)[0:2]
    
video = 1
if video == 1:
    # Set spatial resolution
    x_res = 1000
    y_res = 1000
    xedges = np.linspace(-np.pi, np.pi, num=x_res+1)
    yedges = np.linspace(-np.pi, np.pi, num=y_res+1)
    # Convert to histogram like data
    H  = np.zeros((len(r), y_res,x_res, 3))
    for i in range(len(r)):
        print(i)
        H[i,:,:,2] = np.histogram2d(theta1[:,i], theta0[:,i], bins=(yedges, xedges))[0]
        H[i,:,:,1] = np.histogram2d(theta1[:,i], theta0[:,i], bins=(yedges, xedges))[0]
        H[i,:,:,0] = np.histogram2d(theta1[:,i], theta0[:,i], bins=(yedges, xedges))[0]
    H[H>0] = 1 ; H = H*255; H=np.flip(H,1)
    H=H.astype(dtype='uint8')
    imageio.mimwrite('tipping_circle_map_1.mp4', H , fps = 4)
    

#-------- Degree bifurcation plot beta=4, b=4 ---------------------------------
eps = 0.4
th01, th11, deg1, len1, tip1, F_dot1, lift1 = tipping_circle_map(2.55-eps, 1)
th02, th12, deg2, len2, tip2, F_dot2, lift2 = tipping_circle_map(2.55+eps, 1)
plt.figure()
plt.scatter(th01, th11, s=0.01, color='blue')
plt.plot([5],[5], label='$F_{r_c-\epsilon}$', color='blue')
plt.scatter(th02, th12, s=0.01, color='red')
plt.plot([5],[5], label='$F_{r_c+\epsilon}$', color='red')
plt.ylim([-np.pi, np.pi])
plt.xlim([-np.pi, np.pi])
plt.legend()
plt.ylabel('$F_r(\\theta)$'); plt.xlabel('$\\theta$')


#------- Absolutely value of derivative ---------------------------------------
th01, th11, deg1, len1, tip1, F_dot1, lift1 = tipping_circle_map(3, 1)
plt.plot(th01[1:], np.abs(F_dot1)); plt.ylim([0,50])
