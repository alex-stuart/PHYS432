#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 14:04:43 2022

@author: mac
"""

import numpy as np
import matplotlib.pyplot as pl

# Set up the grid, time and grid spacing, and the sound speed squared
Ngrid = 100
Nsteps = 5000
dt = 0.1
dx = 2.0
cs2 = 1  
gamma = 5/3 #adiabatic index

x = np.arange(Ngrid) * dx # grid
f1 = np.ones(Ngrid) # rho
f2 = np.zeros(Ngrid) # rho x u
f3 = np.ones(Ngrid) # 
u = np.zeros(Ngrid+1) # advective velocity (keep the 1st and last element zero)

def advection(f, u, dt, dx):
    # calculating flux terms
    J = np.zeros(len(f)+1) # keeping the first and the last term zero
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) #update

    return f

# Apply initial Gaussian perturbation
Amp, sigma = 0.5, Ngrid/10
#f1 = f1 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)
#f2 = f2 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)
f3 = f3 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

# plotting
pl.ion()
fig, ax = pl.subplots(1,1)

x1, = ax.plot(x, f1, 'ro')

ax.set_xlim([0, dx*Ngrid+1])
ax.set_ylim([0.8, 2.1])

ax.set_xlabel('x')
ax.set_ylabel('Density')

fig.canvas.draw()

for ct in range(Nsteps):
    # advection velocity at the cell interface
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))

    # update density, momentum and energy
    f1 = advection(f1, u, dt, dx)
    f2 = advection(f2, u, dt, dx)
    f3 = advection(f3, u, dt, dx)

    # add the source term
    f2[1:-1] = f2[1:-1] - 0.5 * (dt / dx) * cs2 * (f1[2:] - f1[:-2])
    f3[1:-1] = f3[1:-1] - 0.5 * (dt / dx) * cs2 * (f2[2:] - f2[:-2])
    
    # correct for source term at the boundary (reflective)
    f2[0] = f2[0] - 0.5 * (dt / dx) * cs2 * (f1[1] - f1[0])
    f2[-1] = f2[-1] - 0.5 * (dt / dx) * cs2 * (f1[-1] - f1[-2])
    f3[0] = f3[0] - 0.5 * (dt / dx) * cs2 * (f2[1] - f2[0])
    f3[-1] = f3[-1] - 0.5 * (dt / dx) * cs2 * (f2[-1] - f2[-2])
    
    cs2 = (1-gamma)*(f3/f1 -(f2**2)/(2*f1**2)) # sound speed squared

    # update the plot
    x1.set_ydata(f1)
    fig.canvas.draw()
    pl.pause(0.001)
    

   
"""
Evolving sound wave
isothermal perturbation
@author: Eve J. Lee
March 14th 2022

import numpy as np
import matplotlib.pyplot as pl

# Set up the grid, time and grid spacing, and the sound speed squared
Ngrid = 100
Nsteps = 5000
dt = 0.1
dx = 2.0
cs2 = 1.0 # sound speed squared

x = np.arange(Ngrid) * dx # grid
f1 = np.ones(Ngrid) # rho
f2 = np.zeros(Ngrid) # rho x u
u = np.zeros(Ngrid+1) # advective velocity (keep the 1st and last element zero)

def advection(f, u, dt, dx):
    # calculating flux terms
    J = np.zeros(len(f)+1) # keeping the first and the last term zero
    J[1:-1] = np.where(u[1:-1] > 0, f[:-1] * u[1:-1], f[1:] * u[1:-1])
    f = f - (dt / dx) * (J[1:] - J[:-1]) #update

    return f

# Apply initial Gaussian perturbation
Amp, sigma = 0.5, Ngrid/10
f1 = f1 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)
f2 = f2 + Amp * np.exp(-(x - x.max()/2) ** 2 / sigma ** 2)

# plotting
pl.ion()
fig, ax = pl.subplots(1,1)

x1, = ax.plot(x, f1, 'ro')

ax.set_xlim([0, dx*Ngrid+1])
ax.set_ylim([0.8, 2.1])

ax.set_xlabel('x')
ax.set_ylabel('Density')

fig.canvas.draw()

for ct in range(Nsteps):
    # advection velocity at the cell interface
    u[1:-1] = 0.5 * ((f2[:-1] / f1[:-1]) + (f2[1:] / f1[1:]))

    # update density and momentum
    f1 = advection(f1, u, dt, dx)
    f2 = advection(f2, u, dt, dx)

    # add the source term
    f2[1:-1] = f2[1:-1] - 0.5 * (dt / dx) * cs2 * (f1[2:] - f1[:-2])

    # correct for source term at the boundary (reflective)
    f2[0] = f2[0] - 0.5 * (dt / dx) * cs2 * (f1[1] - f1[0])
    f2[-1] = f2[-1] - 0.5 * (dt / dx) * cs2 * (f1[-1] - f1[-2])

    # update the plot
    x1.set_ydata(f1)
    fig.canvas.draw()
    pl.pause(0.001)    """