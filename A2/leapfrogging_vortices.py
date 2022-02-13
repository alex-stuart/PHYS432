#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 22:10:09 2022

@author: Alexandre Stuart. Collaborators: Ruijia Yang and Yacine Benkirane
"""

import numpy as np
import matplotlib.pyplot as plt
import imageio

dt = 0.1 #insert appropriate time step 
Nsteps = 130 #insert appropriate total number of timesteps

## Setting up initial conditions (vortex centres and circulation)
# Vortex rings
y_v = [-1,1,-1,1]
x_v = [1,1,3,3]
k_v = [-1,1,-1,1]
#Setting up the plot 
plt.ion()
fig,ax = plt.subplots(1,1)
#mark the initial positions of vortices
p, = ax.plot(x_v,y_v,'*k',markersize=10)
#play around with the marker size and type as you see fit

# draw the initial velocity streamline
ngrid = 2 #dimension of your simulation grid
#Y, X = np.mgrid[-ngrid:ngrid:360j, -ngrid:ngrid:360j]
Y, X = np.mgrid[-4:4:500j, 0:44:500j]
#360j sets the resolution of the cartesian grid; play around with it as you see fit
vel_x = np.zeros(np.shape(X))#this holds x−velocity
vel_y = np.zeros(np.shape(Y))#this holds y−velocity

for i in range(len(x_v)): #looping over each vortex
    #insert lines for computing the total velocity field
    for j in range(len(vel_x[0])):#looping through the field columns
        for k in range(len(vel_x)):#looping through the field rows
            diff_y = y_v[i] - Y[k][j]
            diff_x = x_v[i] - X[k][j]
            #and adding the advection velocity for each point on the field
            if (y_v[i]-Y[k][j])**2+(x_v[i]-X[k][j])**2 !=0:#excluding the advection of vortices on themselves to avoid division by 0
                vel_x[k][j]+=k_v[i]*diff_y/(diff_x**2+diff_y**2)
                vel_y[k][j]+=-k_v[i]*diff_x/(diff_x**2+diff_y**2)       

#setting up the boundaries of the simulation box
ax.set_xlim([0, 44])
ax.set_ylim([-ngrid, ngrid])    

#initial plot of the streamlines
plt.streamplot(X,Y,vel_x,vel_y,density=[1,1],color='b')
plt.savefig('anim_frame0.png')

#Evolution
count=0
while count<Nsteps:
    #re-initializing the vortex advection velocities
    vel_x_vortex=[0,0,0,0]
    vel_y_vortex=[0,0,0,0]
    #and updating them
    for i in range(len(x_v)):
        for j in range(len(x_v)):
            diff_y = y_v[j] - y_v[i] 
            diff_x = x_v[j] - x_v[i]
            if diff_y**2+diff_x**2 !=0:
                vel_x_vortex[i]+=k_v[j]*diff_y/(diff_x**2+diff_y**2)
                vel_y_vortex[i]+=-k_v[j]*diff_x/(diff_x**2+diff_y**2) 
    #updating the positions of vortices
    for i in range(len(x_v)):    
        x_v[i]=x_v[i]+vel_x_vortex[i]*dt
        y_v[i]=y_v[i]+vel_y_vortex[i]*dt

    #re-initializing the total velocity field
    vel_x = np.zeros(np.shape(X))
    vel_y = np.zeros(np.shape(Y))
    #and updating the streamlines 
    for i in range(len(x_v)):
        for j in range(len(vel_x[0])):
            for k in range(len(vel_x)):
                diff_y = y_v[i] - Y[k][j]
                diff_x = x_v[i] - X[k][j]
                if (y_v[i]-Y[k][j])**2+(x_v[i]-X[k][j])**2 !=0:
                    vel_x[k][j]+=k_v[i]*diff_y/(diff_x**2+diff_y**2)
                    vel_y[k][j]+=-k_v[i]*diff_x/(diff_x**2+diff_y**2)
    
    ##update plot
    #the following two lines clear out the previous streamlines 
    ax.collections = []
    ax.patches = []
    
    p.set_xdata(x_v)
    p.set_ydata(y_v)
    
    plot=ax.streamplot(X,Y,vel_x,vel_y,density=[1,1],color='b')
    plt.savefig('anim_frame'+ str(count+1) +'.png')
    count+=1

# Building GIF and saving it as "leapfrog.gif"-code taken from Carvalho, Thiago. Basics of GIFs with Python’s Matplotlib. Available at https://towardsdatascience.com/basics-of-gifs-with-pythons-matplotlib-54dd544b6f30 
    with imageio.get_writer('leapfrog.gif', mode='I') as writer: 
        for filename in ['anim_frame'+ str(count) +'.png' for count in range(121)]:
            image = imageio.imread(filename)
            writer.append_data(image)
      

