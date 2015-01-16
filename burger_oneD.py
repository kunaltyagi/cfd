# -*- coding: utf-8 -*-
"""
Created on Tue Jan  6 13:06:36 2015

@author: kunaltyagi
"""

#!/bin/usr/python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 13:55:28 2014

@author: kunaltyagi
"""

import numpy as np       # NumPy (multidimensional arrays, linear algebra, ...)
import scipy as sp       # SciPy (signal and image processing library)
import sympy
import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import *              # Matplotlib's pylab interface
ion()                            # Turned on Matplotlib's interactive mode

x, nu, t = sympy.symbols('x nu t')
phi = sympy.exp(-(x-4*t)**2/(4*nu*(t+1))) + sympy.exp(-(x-4*t-2*np.pi)**2/(4*nu*(t+1)))
phiprime = phi.diff(x)
print phiprime
u = -2*nu*(phiprime/phi)+4
ufunc = sympy.utilities.lambdify((t, x, nu), u)

nx = 101
nt = 300
dx = 2*np.pi/(nx-1)
nu = .07
dt = dx*nu

x = np.linspace(0, 2*np.pi, nx)
#u = np.empty(nx)
un = np.empty(nx)
t = 0

#u = np.asarray([ufunc(t, x0, nu) for x0 in x])

u = np.asarray([sin(x0) for x0 in x])

#plot(u, 'r')
u0 = np.asarray([sin(x0) for x0 in x])

for n in range(nt):
    un = u0.copy()
    for i in range(nx-1):
        if u0[i] > 0:
            u0[i] = un[i] - dt/dx *(un[i]**2 - un[i-1]**2)/2
        else:
            u0[i] = un[i] - dt/dx *(un[i+1]**2 - un[i]**2)/2
            
    u0[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2])

#plot(u)

u1 = np.asarray([sin(x0) for x0 in x])

for n in range(nt):
    un = u1.copy()
    for i in range(nx-1):
        if u1[i] > 0:
            u1[i] = un[i] - un[i]*dt/dx *(un[i] - un[i-1])
        else:
            u1[i] = un[i] - un[i]*dt/dx *(un[i+1] - un[i])
    u1[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2])

#plot(u1, 'y')
"""        
u_analytical = np.asarray([ufunc(nt*dt, xi, nu) for xi in x])
"""
u2 = np.asarray([sin(x0) for x0 in x])
for n in range(nt):
    un = u2.copy()
    for i in range(nx-1):
        if u2[i] > 0:
            u2[i] = un[i] - un[i] * dt/dx *(un[i] - un[i-1]) + nu*dt/dx**2*\
                (un[i+1]-2*un[i]+un[i-1])
        else:
            u2[i] = un[i] - un[i] * dt/dx *(un[i+1] - un[i]) + nu*dt/dx**2*\
                (un[i+1]-2*un[i]+un[i-1])            
    u2[-1] = un[-1] - un[-1] * dt/dx * (un[-1] - un[-2]) + nu*dt/dx**2*\
                (un[0]-2*un[-1]+un[-2])
plot(u, 'b', u0, 'r', u1, 'y', u2,'g')

legend(('original', "w/o nu, u*(u')", "w/o nu, ((u**2)/2)'", 'viscous'))