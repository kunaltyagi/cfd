#!/bin/usr/python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 13:55:28 2014

@author: kunaltyagi
"""

import numpy as np       # NumPy (multidimensional arrays, linear algebra, ...)
# import scipy as sp       # SciPy (signal and image processing library)
import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import *              # Matplotlib's pylab interface
ion()                            # Turned on Matplotlib's interactive mode

# @systemProperties
TIME = 0.00075
DELTA_X = 0.01
CFL = 0.9
R = 287.0
GAMMA = 1.4

# @derivedConstants
NUM_X = (int)(1/DELTA_X + 2)

# @equation
# \partialDerivative{U}{t} + \partialDerivative{F}{x} = 0

# @variables
U       = np.array(0.0)
F       = np.array(0.0)
F_plus  = np.array(0.0)
F_minus = np.array(0.0)
u   = np.array(0.0)
M   = np.array(0.0)
rho = np.array(0.0)
p   = np.array(0.0)
T   = np.array(0.0)
a   = np.array(0.0)
SR  = 0.0   # spectral radius
P0  = 101325.0

U.resize([NUM_X, 3])
F.resize([NUM_X, 3])
F_plus.resize([NUM_X, 3])
F_minus.resize([NUM_X, 3])
u.resize([NUM_X])
M.resize([NUM_X])
rho.resize([NUM_X])
p.resize([NUM_X])
T.resize([NUM_X])
a.resize([NUM_X])

# @supplememtry variables
F1 = zeros((NUM_X,3))
F2 = zeros((NUM_X,3))
F3 = zeros((NUM_X,3))

# @initialConditions
initial_temp = 300
T = initial_temp * np.ones([NUM_X])
for x in range(0, NUM_X):
    p[x]   = 5.0*P0 if x < (0.5 * NUM_X) else P0
    rho[x] = p[x]/(R * T[x])
    U[x][:] = [rho[x], rho[x]*u[x], rho[x]*(u[x]**2)/2.0 + p[x]/(GAMMA-1)]
    
# @algorithm
t = 0
while t < TIME:
    SR = 0
    for x in range(1, NUM_X-1):
        a[x] = sqrt(GAMMA*R*T[x])
        
        # @fluxVectorSplitting        
        velPar = (u[x]**2)/2.0 + (a[x]**2)/(GAMMA-1)
        d = rho[x]*(u[x]+a[x])/(2*GAMMA)
        e = (GAMMA-1)*rho[x]*u[x]/GAMMA
        f = rho[x]*(u[x] - a[x])/(2*GAMMA)
        
        F1[x][:] = [d, d*(u[x]+a[x]), d*(velPar + a[x]*u[x])]
        F2[x][:] = [e, e*u[x],        e*u[x]**2/2.0]
        F3[x][:] = [f, f*(u[x]-a[x]), f*(velPar - a[x]*u[x])]
       
        M[x] = u[x]/a[x]
        if M[x] < -1:
            F_plus[x][:]  = zeros(3)
            F_minus[x][:] = F1[x][:] + F2[x][:] + F3[x][:]
        elif M[x] <= 0:
            F_plus[x][:]  = F1[x][:]
            F_minus[x][:] = F2[x][:] + F3[x][:]
        elif M[x] <=1:
            F_plus[x][:]  = F2[x][:] + F1[x][:]
            F_minus[x][:] = F3[x][:]
        elif M[x] > 1:
            F_plus[x][:]  = F1[x][:] + F2[x][:] + F3[x][:]
            F_minus[x][:] = zeros(3)
            
        # @timeStep for stability condition
        temp = max(abs(u[x]),abs(u[x]+a[x]),abs(u[x]-a[x]))
        if temp > SR:
            SR = temp
            
    # @leftBoundaryCondition
    F_plus[0][:] = F_plus[1][:]
    F_minus[0][:] = F_minus[1][:]

    # @rightBoundaryCondition
    F_plus[-1][:] = F_plus[-2][:]
    F_minus[-1][:] = F_minus[-2][:]
    
    # @nextIteration
    delta_t = (CFL*DELTA_X)/SR
    alpha = delta_t/DELTA_X
    t = t + delta_t
    

    # @calculating U
    for x in range(1, NUM_X-1):
        U[x][:] = U[x][:] - alpha*(F_plus[x][:] - F_plus[x-1][:]) - alpha*(
                  F_minus[x+1][:] - F_minus[x][:])
        rho[x] = U[x][0]
        u[x] = U[x][1]/rho[x]
        p[x] = (GAMMA-1)*(U[x][2] - (U[x][1]**2)/(2.0*U[x][0]))
        T[x] = p[x]/(R*rho[x])
        a[x] = sqrt(GAMMA*R*T[x])
        M[x] = u[x]/a[x]
        
u_max = max(u)
rho_max = max(rho)
p_max = max(p)

# Normalized velocity
plot(u/u_max,'r')
# Normalized density
plot(rho/rho_max, 'b')
# Normalized Pressure
plot(p/p_max, 'g')
# Mach number
plot(M, 'y')