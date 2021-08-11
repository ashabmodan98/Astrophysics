# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 11:54:09 2019

@author: DELL
"""


import matplotlib.pyplot as plt
from math import cos, sin
import numpy as np
# Constants
M_sun = 1.9891 * 10**30
G = 6.67408 * 10** -11
PI = np.pi

m1 = float(input("Enter the mass of star 1(in solar mass):"))
m2 = float(input("Enter the mass of star 2(in solar mass):"))
e = float(input("Enter the orbital eccentricity (<1):"))
P = float(input("Enter the orbital timeperiod (in days):"))

# Unit Conversion

m1 *= M_sun
m2 *= M_sun
P *= 24 * 3600

# Reduced mass system
M = m1 + m2
m = m1 * m2 / M  
a = (G * M * (P ** 2) / (4 * PI**2))**(1 / 3)


# Semi-Major axes
a1 = m * a / m1
a2 = m * a / m2


# Inital Conditions
t = []	   #time
v1r = []
v2r = []

t.append(0)
dt = P / 10000 # Rate of change of time
theta = 0
L = m * (G * M * a * (1 - e**2))**0.5 # Angular Momentum
dA = L / (2 * m)

i = 0
while t[i] + dt < 3.5 * P:              # Time for 1 revolution
    r = a * (1 - e**2) / (1 + e * cos(theta))		# position of CM
    vr = -sin(theta) * (G * M * (2 / r - 1 / a))**0.5  # Radial velocity of CM

    # Required values of radial velocity, Y axis
    v1r.append(vr * m / (m1 * 1000))
    v2r.append(-vr * m / (m2 * 1000))

    dtheta = 2 * dA * dt * 1 / (r**2)
    theta += dtheta
    t.append(t[i] + dt)
    i += 1
t = [i / (24 * 3600) for i in t]   # in Days

plt.plot(t[:-1], v1r, 'b')
plt.plot(t[:-1], v2r, 'r')
plt.xlabel('Time (Days)')
plt.ylabel('Radial Velocity  (km/sec)')
plt.legend(('Star 1','Star 2'),loc='upper right')
plt.show()