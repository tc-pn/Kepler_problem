"""
2D Kepler problem solver (i.e. 2 objets orbiting each other).
I intend to solve for any masses.
"""
import numpy as np
import pandas as pd
import pickle

from KeplerObjectClass import KeplerObject

def main(Object1 : KeplerObject,
         Object2 : KeplerObject,
         nsteps : int=100) -> pd.DataFrame:
    """
    Main function of the solver.
    It is a standard RK2 procedure.
    Parameters:
    Object1 : KeplerObject
        First object in the system
    Object2 : KeplerObject
        Second object in the system
    nsteps: int
        Number of time steps over which to solve the equations of motion
    """
    delta_t = 1./365.24
    Object1, Object2 = RK2_procedure(Object1=Object1,
                               Object2=Object2,
                               nsteps=nsteps,
                               delta_t=delta_t)
    return



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

    Sun = KeplerObject(mass=m1,
                       x_pos=x1,
                       y_pos=x2,
                       vx_ini=vx1_ini,
                       vy_ini=vy1_ini)
    Earth = KeplerObject(mass=m2,
                         x_pos=x2,
                         y_pos=y2,
                         vx_ini=vx2_ini,
                         vy_ini=vy2_ini)



