#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 22:01:34 2020

@author: robertdenomme
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from tqdm import tqdm
from WaveEngine import WaveEngine
import cmath


# fig = plt.figure()
# ax = fig.add_subplot(121, projection='3d',
#                       zlim=(-1, 1))
# ax2 = fig.add_subplot(122, projection='3d',
#                        zlim=(-1, 1))


fig, (ax,ax2) = plt.subplots(1,2,subplot_kw=dict(projection='3d', zlim=(-1,1)))
#can also use a instead of (ax,ax2) to get multiple axis return objects


Nfrm =10*20
fps = 10
alpha=2

x_samples = 40
y_samples = x_samples

xs = np.linspace(-1, 1, x_samples)
delta_x = 2/x_samples
ys = np.linspace(-1, 1, y_samples)
delta_y = 2/y_samples

X, Y = np.meshgrid(xs, ys)

def init_bump(X,Y,x_0=0, y_0=0):
    radius=20
    height=.5
    W=0*X
    W=height*np.exp(-radius*((X-x_0)**2+(Y-y_0)**2))
    return W


def plot_from_initial_state(V, V_dot, Wave, axi):
    """Get the initial conditions of position and momentum V, V_dot and 
    plot the evolution of the system."""
    Z,Z_dot = ([V],[V_dot])
    ims=[]
    for i in tqdm(range(Nfrm)):
        W,W_dot = Wave.get_next(Z[-1],Z_dot[-1])
        Z.append(W)
        Z_dot.append(W_dot)
    for i in tqdm(range(Nfrm)):
        im = axi.plot_wireframe(X, Y, Z[i], rstride=1, cstride=1, color='k', 
                                linewidth=0.5)
        ims.append([im])
    return ims
        
            
def get_mode_initial_state(Lambda, E, i):
    """Get two vectors that form an oscillatory state. upon return, 
    the evolution of the system is cos(omega*t)Re+sin(omega*t)Im"""
    Re,_ = Wave.from_phase_space(E[i].real)
    Im,_ = Wave.from_phase_space(E[i].imag)
    #note, we don't need the dot since these belong to the same 
    #complex eigenvalue
    omega = cmath.phase(Lambda[i])
    return (Re,Im,omega)

def plot_mode(Re,Im,omega,axi):
    """Use the eigenvalue propery to plot it wouthout the dynamics 
    matrix explicitly"""
    Z=[]
    ims = []
    
    for i in range(Nfrm):
        U=np.cos(omega*i)*Re+np.sin(omega*i)*Im
        Z.append(U)
    for i in tqdm(range(Nfrm)):
        
        im = axi.plot_wireframe(X, Y, Z[i], rstride=1, cstride=1, color='k', 
                                linewidth=0.5)
        
        # frame_text = ax.text(-1, -1, 1, '', transform=ax.transAxes)
        # frame_text.set_text(f"frame: {i}")

        ims.append([im])
    return ims

def plot_from_trajectory_data(Z, axi):
    ims = []
    
    for i in range(Nfrm):
        im = axi.plot_wireframe(X, Y, Z[i], rstride=1, cstride=1, color='k', 
                                linewidth=0.5)
        ims.append([im])
    
    return ims


M,N = np.shape(X)
Wave = WaveEngine(np.shape(X),Nfrm, alpha=alpha)
Z,Z_dot=([],[])
# Z.append(init_bump(X,Y))
# Z_dot.append(np.zeros(np.shape(X)))
i=1
Lambda, E = np.linalg.eig(Wave.L)
V = Wave.unflatten(E[:,99].real)
V_dot = 0*V
ims1 = plot_from_initial_state(V, V_dot, Wave, ax)


# W = init_bump(X,Y)
# W_dot = 0*W
# ims2 = plot_from_initial_state(W, W_dot, Wave, ax2)


W = init_bump(X,Y)
W_dot = 0*W
Z = Wave.get_continuous_trajectory(W,W_dot)
ims2 = plot_from_trajectory_data(Z,ax2)

ims=[im1+im2 for im1,im2 in zip(ims1,ims2)]

"""
ones_vector = np.zeros((M*N,1),dtype=float)+1
ones_vector.transpose()
Signal = []
for X in Z:
    Signal.append(ones_vector@X@ones_vector)
T = [i for i in range(Signal)]

figgy, axy = plt.subplots()
axy.plot(T,Signal)
# ims=ims1


"""


ani = animation.ArtistAnimation(fig, ims,interval=1000/fps, blit=True, repeat_delay=0)

print("animations finished")

# fn = 'plot_wireframe_funcanimation'
# ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)

print("files saved")

