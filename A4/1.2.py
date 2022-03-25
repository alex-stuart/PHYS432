#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 10:05:26 2022

@author: mac
"""

import matplotlib.pyplot as plt
import numpy as np

m1=0.01
m2=0.1
m3=1
m4=10
m5=100

U = np.logspace(-1, 3, 10000)
e1=2.3*10**(-2)*(m1/0.035)**(-1/3)*(U/10)**2+2.54*10**(-2)*(m1/0.035)**(1/3)*(U/10)**(-2)+1.077*(m1/0.035)**(-1/4)*(U/10)**(-1)
e2=2.3*10**(-2)*(m2/0.035)**(-1/3)*(U/10)**2+2.54*10**(-2)*(m2/0.035)**(1/3)*(U/10)**(-2)+1.077*(m2/0.035)**(-1/4)*(U/10)**(-1)
e3=2.3*10**(-2)*(m3/0.035)**(-1/3)*(U/10)**2+2.54*10**(-2)*(m3/0.035)**(1/3)*(U/10)**(-2)+1.077*(m3/0.035)**(-1/4)*(U/10)**(-1)
e4=2.3*10**(-2)*(m4/0.035)**(-1/3)*(U/10)**2+2.54*10**(-2)*(m4/0.035)**(1/3)*(U/10)**(-2)+1.077*(m4/0.035)**(-1/4)*(U/10)**(-1)
e5=2.3*10**(-2)*(m5/0.035)**(-1/3)*(U/10)**2+2.54*10**(-2)*(m5/0.035)**(1/3)*(U/10)**(-2)+1.077*(m5/0.035)**(-1/4)*(U/10)**(-1)

plt.xscale('log')
plt.yscale('log')
plt.plot(U,e1)
plt.plot(U,e2)
plt.plot(U,e3)
plt.plot(U,e4)
plt.plot(U,e5)
plt.title('Cost of Transport e_cost as a function of flight speed for five birds of different masses')
plt.legend(['m=0.01kg', 'm=0.1kg', 'm=1kg', 'm=10kg', 'm=100kg'])
plt.xlabel('Bird flight speed (m/s)')
plt.ylabel('Cost of Transport e_cost (unitless)')
plt.show()

#find minima
min_U1 = U[e1.argmin()]
min_U2 = U[e2.argmin()]
min_U3 = U[e3.argmin()]
min_U4 = U[e4.argmin()]
min_U5 = U[e5.argmin()]

print(min_U1, min(e1))
print(min_U2, min(e2))
print(min_U3, min(e3))
print(min_U4, min(e4))
print(min_U5, min(e5))
