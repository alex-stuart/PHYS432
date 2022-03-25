#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:55:33 2022

@author: mac
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


m=np.array([0.01,0.1,1,10,100])
e_min=np.array([0.802,0.426,0.230,0.131,0.083])

def fit_func(m, a, b):
    return a*m**b

params = curve_fit(fit_func, m, e_min)

[a, b] = params[0]
print([a, b])

mvals = np.logspace(-2, 2, 3000)

plt.scatter(m,e_min,marker='o',color='r')
plt.plot(mvals, fit_func(mvals, a, b))
plt.xscale('log')
plt.yscale('log')
plt.title('Minimum Cost of Transport e_cost as a function of bird mass')
#plt.legend(['m=0.01kg', 'm=0.1kg', 'm=1kg', 'm=10kg', 'm=100kg'])
plt.xlabel('Bird mass (kg)')
plt.ylabel('Minimum cost of Transport e_cost (unitless)')
plt.show()