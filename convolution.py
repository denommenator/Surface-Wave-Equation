#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 12:38:11 2020

@author: robertdenomme
"""

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
print('numpy: '+np.version.full_version)
import matplotlib.animation as animation
import matplotlib
print('matplotlib: '+matplotlib.__version__)
from matplotlib.colors import LightSource 
from matplotlib import cm
from scipy import ndimage
#import scipy.ndimage.filters


Nfrm = 24
fps = 24
alpha=2


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Make the X, Y meshgrid.
xs = np.linspace(-1, 1, 40)
ys = np.linspace(-1, 1, 40)
X, Y = np.meshgrid(xs, ys)




def laplace(a_frame):
    k = np.array([[0,.25,0],
                  [.25,-1,.25],
                  [0,.25,0]])
    W = ndimage.convolve(a_frame,k,mode='constant', cval=0.0)
    return W
            
def generate_all_frames(Z):
    i=2
    while i<Nfrm:
        Z.append(Z[-1]+(Z[-1]-Z[-2])+alpha*laplace(Z[-1]))
        i+=1
    return Z

def init_bump(X,Y,x_0=0, y_0=0):
    radius=20
    height=.5
    W=0*X
    W=height*np.exp(-radius*((X-x_0)**2+(Y-y_0)**2))
    return W
            
      
Z = []
Z.append(init_bump(X,Y))
Z.append(Z[0])
# print(Z[0])
# print(Z[1])
#print (Z) 
generate_all_frames(Z)
#print (Z[2])
print("frames generated")


# Set the z axis limits so they aren't recalculated each frame.
ax.set_zlim(-1, 1)

wframe = None

def update(idx):
    global wframe
    # If a line collection is already remove it before drawing.
    if wframe:
        ax.collections.remove(wframe)

    # Plot the new wireframe and pause briefly before continuing.
    W = np.array(Z[idx])
    #print(f"Frame: {idx} shift: {.05*idx}")
    ls = LightSource(45,45)
    rgb = np.ones((W.shape[0], W.shape[1], 3))
    illuminated_surface = ls.shade_rgb(rgb,W)
    wframe = ax.plot_surface(X, Y, W, rstride=1, cstride=1, color='k', linewidth=0, facecolors=illuminated_surface)

ani = animation.FuncAnimation(fig, update, Nfrm, interval=1000/fps)




# fn = 'plot_wireframe_funcanimation'
# ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)

