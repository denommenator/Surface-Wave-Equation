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

Nfrm = 400
fps = 30
alpha=1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Make the X, Y meshgrid.
xs = np.linspace(-1, 1, 50)
ys = np.linspace(-1, 1, 50)
X, Y = np.meshgrid(xs, ys)



def laplace(a_frame):
    W=[]
    for k in range(len(a_frame)):
        W.append([])
        for l in range(len(a_frame[0])):
            if k==0:
                n=0
            else:
                n=a_frame[k-1][l]
            if k==len(a_frame)-1:
                s=0
            else:
                s=a_frame[k+1][l]
            if l==0:
                w=0
            else:
                w=a_frame[k][l-1]
            if l<len(a_frame[0])-1:
                e=a_frame[k][l+1]
            else:
                e=0
            W[k].append(1/4*(n+s+e+w) - a_frame[k][l])
    return W
            
def generate_all_frames(Z):
    i=2
    while i<Nfrm:
        Z.append(Z[-1]+Z[-1]-Z[-2]+alpha*laplace(Z[-1]))
        i+=1
    return Z

           
      
Z = []
Z.append(X*0+1/4)
Z.append(Z[0])
#print (Z) 
generate_all_frames(Z)
#print (Z[2])



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
    wframe = ax.plot_wireframe(X, Y, W, rstride=1, cstride=1, color='k', linewidth=0.1)

ani = animation.FuncAnimation(fig, update, Nfrm, interval=1000/fps)




fn = 'plot_wireframe_funcanimation'
ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)

