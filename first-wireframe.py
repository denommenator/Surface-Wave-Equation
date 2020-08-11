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

Nfrm = 10
fps = 10

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Make the X, Y meshgrid.
xs = np.linspace(-1, 1, 100)
ys = np.linspace(-1, 1, 100)
X, Y = np.meshgrid(xs, ys)

Z = np.sin(10*(X**2+Y**2))




# Set the z axis limits so they aren't recalculated each frame.
ax.set_zlim(-1, 1)

wframe = None

def update(idx):
    global wframe, frame_no
    # If a line collection is already remove it before drawing.
    if wframe:
        ax.collections.remove(wframe)

    # Plot the new wireframe and pause briefly before continuing.
    W = Z+idx*.05
    print(f"Frame: {idx} shift: {.05*idx}")
    wframe = ax.plot_wireframe(X, Y, W, rstride=1, cstride=1, color='k', linewidth=0.1)

ani = animation.FuncAnimation(fig, update, Nfrm, interval=1000/fps)




# fn = 'plot_wireframe_funcanimation'
# ani.save(fn+'.mp4',writer='ffmpeg',fps=fps)

