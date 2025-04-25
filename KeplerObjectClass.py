"""
Module defining the KeplerObject class.
"""
import numpy as np


class KeplerObject():
    """
    Class describing the mass and positions on the (x,y) plane of a massive point-like object.
    """
    def __init__(self,
                 mass : float,
                 x_pos : float,
                 y_pos : float,
                 vx : float,
                 vy : float):
        self.G = 39.478 # AU**3 * year**-2. / solar mass
        self.mass = mass # solar mass
        self.x_pos = x_pos # AU
        self.y_pos = y_pos # AU
        self.vx = vx # AU / yr
        self.vy = vy # AU / yr
        self.positions = [[self.x_pos, self.y_pos]]
        self.velocities = [[self.vx, self.vy]]

    def update_vector(self,
                      new_x : float,
                      new_y : float,
                      new_vx : float,
                      new_vy : float):
        self.x_pos = new_x
        self.y_pos = new_y
        self.vx = new_vx
        self.vy = new_vy
        self.positions.append([self.x_pos, self.y_pos])
        self.velocities.append([self.vx, self.vy])

    def acceleration(self,
                     rel_pos : float,
                     mass_other : float,
                     x_pos_other : float,
                     y_pos_other : float) -> float:
        d_x = np.abs(self.x_pos - x_pos_other)
        d_y = np.abs(self.y_pos - y_pos_other)
        out = - self.G * mass_other * rel_pos / (d_x**2. + d_y**2.)**(3./2.)
        return out