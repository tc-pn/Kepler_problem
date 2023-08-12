# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 19:46:29 2023

@author: Thomas Czuba
"""

import numpy as np
import matplotlib.pyplot as plt

G = 4. * np.pi**2. # in astronomical units, years and solar mass
dt = 10e-5 # in years 


masses = np.array([1.,1.]) # in solar masses

#sun then earth
r = np.array([1.,1.]) # in astronomical units
theta = np.array([np.pi,0.])

v_r = np.array([0.,0.]) # in astronomical units per year
#v_theta = np.array([0., 2. * np.pi]) # in radians per year
v_theta = np.array([1., 1.]) * np.sqrt(G * np.sum(masses)) / r**(3./2.) 

J = r**2. * v_theta

class Kepler:
    def __init__(self,r,theta,v_r,v_theta,masses):
        self.r = r
        self.theta = theta
        self.v_r = v_r
        self.v_theta = v_theta
        self.masses = masses
        
    def potential_without_mass(self, positions, theta):
        norm = (positions[0]**2. + positions[1]**2. 
                       - 2. * positions[1] * positions[0] * np.cos(np.abs(theta[1]-theta[0])))
        return - G / norm**2. + J**2. / norm**3.
        
    def evolution(self):
        v_r_foo = self.v_r
        r_foo = self.r
        
        v_r_foo = v_r_foo + dt / 2. * self.masses * self.potential_without_mass(r,theta)
        r_foo = r_foo + dt / 2. * v_r_foo
        
        v_theta_foo = J / r_foo**2.
        theta_foo = self.theta + dt / 2. * v_theta_foo
        
        self.v_r = self.v_r + dt * self.masses * self.potential_without_mass(r_foo,theta_foo)
        self.r = self.r + dt * v_r_foo
        
        self.v_theta = J / self.r**2.
        self.theta = self.theta + dt * self.v_theta
        
        pass
    pass

Test = Kepler(r,theta,v_r,v_theta,masses)
positions = []
angles = []
for i in range(21400):
    Test.evolution()
    positions.append(Test.r)
    angles.append(Test.theta)

positions = np.array(positions)
angles = np.array(angles)

#print(positions)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(angles[:,0],positions[:,0], label = 'Sun', color = 'blue', linewidth = 1.)
ax.plot(angles[:,1],positions[:,1], label = 'Earth', color = 'red', linewidth = 1.)
ax.scatter(angles[-1,0],positions[-1,0], color = 'blue')
ax.scatter(angles[-1,1],positions[-1,1], color = 'red')

#fig, ax = plt.subplots()

#time = np.arange(0,30000) * dt

#ax.plot(time,positions[:,0], label = 'Sun', color = 'blue', linewidth = 2.)
#ax.plot(time,angles[:,0], label = 'Sun', color = 'blue', linewidth = 2.)
#ax.plot(time,positions[:,1], label = 'Earth', color = 'red', linewidth = '2.')
#ax.plot(time,angles[:,1], label = 'Earth', color = 'red', linewidth = '2.')
plt.show()
