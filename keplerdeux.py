# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 19:46:29 2023

@author: Thomas Czuba
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

G = 4. * np.pi**2. # in astronomical units, years and solar mass
dt = 1. / 365. * 0.01

masses = np.array([3. * 10.**(-6.),1.]) # in solar masses

#sun then earth
r = np.array([0.00000001,1.]) # in astronomical units
theta = np.array([0.,0.])
v_r = np.array([0.,0.])
v_theta = np.array([0., 1.]) * np.sqrt(G * np.sum(masses)) / r**(3./2.) 


J = r**2. * v_theta

class Kepler:
    def __init__(self,r,theta,v_r,v_theta,masses,Ntime):
        self.r = r
        self.theta = theta
        self.v_r = v_r
        self.v_theta = v_theta
        self.masses = masses
        self.Ntime = Ntime
        
    def force(self, positions, angle, mass):
        diff_angles = self.theta[1] - self.theta[0]
        norm = np.sqrt(positions[0]**2. + positions[1]**2. 
                       - 2. * positions[1] * positions[0] * np.cos(diff_angles))
        return - G * mass / norm**2. + J**2. / norm**3.
    
    def energy(self):
        diff_angles = self.theta[1] - self.theta[0]
        norm = np.sqrt(self.r[0]**2. + self.r[1]**2. 
                       - 2. * self.r[1] * self.r[0] * np.cos(diff_angles))
        kin_r = 0.5 * (self.masses[1] * self.v_r[0]**2. + self.masses[0] * v_r[1]**2.)
        
        pot = - G * self.masses[0] * self.masses[1] / norm + np.sum(J**2. / (2. * norm**3.))
        return kin_r + pot
    
    def evolution_RKTWO(self):
        
        J = self.r**2. * self.v_theta

        v_r_foo = self.v_r
        r_foo = self.r
        theta_foo = self.theta
        
        v_r_foo = v_r_foo + dt / 2. * self.force(r_foo, theta_foo, masses)
        r_foo = r_foo + dt / 2. * v_r_foo
        
        v_theta_foo = J / r_foo**2.
        theta_foo = self.theta + dt / 2. * v_theta_foo
        
        self.v_r = self.v_r + dt * self.force(r_foo, theta_foo, masses)
        self.r = self.r + dt * v_r_foo
        
        self.v_theta = J / self.r**2.
        self.theta = self.theta + dt * self.v_theta

        pass
    
    def simp_integrator(self):
        
        J = self.r**2. * self.v_theta
        
        acc_t = self.force(self.r, self.theta, self.masses)
        
        v_r_foo = self.v_r + dt / 2. * acc_t
        self.r = self.r + v_r_foo * dt
        
        self.v_theta = J / self.r**2.
        self.theta = self.theta + dt * self.v_theta
        
        acc_t_plus_dt = self.force(self.r, self.theta, self.masses)
        
        self.v_r = v_r_foo + dt / 2. * acc_t_plus_dt
        
        pass
    
    
    def solver(self):
        E = []
        radial_pos, angles, velocity_r, velocity_theta = [], [], [], []
        for i in range(self.Ntime):
            self.simp_integrator()
            if i % 100 == 0:
                print(dt * i, dt, i)
                radial_pos.append(self.r)
                angles.append(self.theta)
                E.append(self.energy())
                velocity_r.append(self.v_r)
                velocity_theta.append(self.v_theta)
                
        radial_pos, angles = np.array(radial_pos), np.array(angles)
        velocity_r, velocity_theta = np.array(velocity_r), np.array(velocity_theta)

        
        df_radial = pd.DataFrame(radial_pos, columns = ["r_one", "r_two"])
        df_angles = pd.DataFrame(angles, columns = ["angle_one", "angle_two"])
        df_energy = pd.DataFrame({"Time": np.arange(0,Ntime,100) * dt, "Energy": E})
        df_velocity_r = pd.DataFrame(velocity_r, columns = ["v_r_one", "v_r_two"])
        df_velocity_theta = pd.DataFrame(velocity_theta, columns = ["v_theta_one", "v_theta_two"])

        df_radial.to_csv('radial_distances.csv')
        df_angles.to_csv('angles.csv')
        df_velocity_r.to_csv('radial_velocities.csv')
        df_velocity_theta.to_csv('theta_velocities.csv')
        df_energy.to_csv('energy.csv')
        pass

    def polar_viz(self):
        
        positions = pd.read_csv('radial_distances.csv')
        angles = pd.read_csv('angles.csv')
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(angles["angle_one"],positions["r_one"], label = 'Sun', color = 'blue')
        ax.plot(angles["angle_two"],positions["r_two"], label = 'Earth', color = 'red')
        ax.plot(angles["angle_one"].iloc[-1],positions["r_one"].iloc[-1], color = 'blue', marker = 'o')
        ax.plot(angles["angle_two"].iloc[-1],positions["r_two"].iloc[-1], color = 'red', marker = 'o')
        plt.show()
        pass

    def energy_viz(self):
        ener = pd.read_csv('energy.csv')

        fig, ax = plt.subplots()
        plt.xlabel(r'Time $({\rm yr})$', fontsize=15)
        plt.ylabel(r'E $({\rm M}_\odot {\rm AU}^2 {\rm yr}^{-2})$', fontsize=15)
        ax.plot(ener["Time"],(ener["Energy"] - ener.loc[0,"Energy"]), color = 'blue')
        
        plt.show()
        pass
    
    def pos_vs_time(self):
        positions = pd.read_csv('radial_distances.csv')
        angles = pd.read_csv('angles.csv')
        ener = pd.read_csv('energy.csv')

        nrows, ncolumns = 2, 2
        fig, axs = plt.subplots(nrows, ncolumns, sharex=False, sharey=False)
        axs[0,0].plot(ener["Time"],positions["r_one"], color = 'blue')
        axs[0,1].plot(ener["Time"],positions["r_two"], color = 'blue')
        axs[1,0].plot(ener["Time"],angles["angle_one"], color = 'blue')
        axs[1,1].plot(ener["Time"],angles["angle_two"], color = 'blue')
        plt.show()
        pass

    def v_vs_time(self):
        v_r = pd.read_csv('radial_velocities.csv')
        v_theta = pd.read_csv('theta_velocities.csv')
        ener = pd.read_csv('energy.csv')

        nrows, ncolumns = 2, 2
        fig, axs = plt.subplots(nrows, ncolumns, sharex=False, sharey=False)
        axs[0,0].plot(ener["Time"],v_r["v_r_one"], color = 'blue')
        axs[0,1].plot(ener["Time"],v_r["v_r_two"], color = 'blue')
        axs[1,0].plot(ener["Time"],v_theta["v_theta_one"], color = 'blue')
        axs[1,1].plot(ener["Time"],v_theta["v_theta_two"], color = 'blue')
        plt.show()
        pass
    
    def phase_diagrams(self):
        positions = pd.read_csv('radial_distances.csv')
        angles = pd.read_csv('angles.csv')
        v_r = pd.read_csv('radial_velocities.csv')
        v_theta = pd.read_csv('theta_velocities.csv')
        
        nrows, ncolumns = 2, 2
        fig, axs = plt.subplots(nrows, ncolumns, sharex=False, sharey=False)
        for ax in axs.flat:
            ax.set(xlabel='x-label', ylabel='y-label')
        axs[0,0].plot(positions["r_one"],v_r["v_r_one"], color = 'blue')
        axs[0,1].plot(positions["r_two"],v_r["v_r_two"], color = 'blue')
        axs[1,0].plot(angles["angle_one"],v_theta["v_theta_one"], color = 'blue')
        axs[1,1].plot(angles["angle_two"],v_theta["v_theta_two"], color = 'blue')

        plt.show()
        pass

    pass

#Ntime = 10**5
Ntime = 365 * 500
Test = Kepler(r,theta,v_r,v_theta,masses,Ntime)
Test.solver()
Test.pos_vs_time()
Test.polar_viz()
Test.energy_viz()
Test.v_vs_time()
Test.phase_diagrams()