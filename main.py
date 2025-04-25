"""
2D Kepler problem solver (i.e. 2 objets orbiting each other).
I intend to solve for any masses.
"""
import pickle

import numpy as np
import pandas as pd

from SystemKeplerClass import KeplerSystem

if __name__ == "main":
    #G = 6.67e-11
    G = 39.478 # AU**3 * year**-2. / solar mass
    #m_sun = 1.989e30
    m_sun = 1. #solar mass units
    m1 = 1. * m_sun
    m2 = 3.00274e-6 * m_sun
    x1, y1 = 0., 0.
    vx1_ini, vy1_ini = 0., 0.
    #x2, y2 = 150.35e6, 0.
    x2, y2 = 1.0167, 0. # earth aphelion in astronomical units
    vx2_ini, vy2_ini = 0., 6.1744 # earth velocity at aphelion

    System = KeplerSystem(m1=m1, m2=m2,
                          x1_pos=x1, x2_pos=x2,
                          y1_pos=y1, y2_pos=y2,
                          vx1=vx1_ini, vx2=vx2_ini,
                          vy1=vy1_ini, vy2=vy2_ini)




