"""
Module defining the KeplerSystem class.
"""
import numpy as np

from KeplerObjectClass import KeplerObject

class KeplerSystem():
    def __init__(self,
                 m1 : float, m2 : float,
                 x1_pos : float, x2_pos : float,
                 y1_pos : float, y2_pos : float,
                 vx1 : float, vx2 : float,
                 vy1 : float, vy2 : float):
        
        self.G = 39.478 # AU**3 * year**-2. / solar mass
        self.Object1 = KeplerObject(mass=m1,
                                    x_pos=x1_pos,
                                    y_pos=y1_pos,
                                    vx=vx1,
                                    vy=vy1)
        self.Object2 = KeplerObject(mass=m2,
                                    x_pos=x2_pos,
                                    y_pos=y2_pos,
                                    vx=vx2,
                                    vy=vy2)
        self.x1_ini = x1_pos #AU
        self.x2_ini = x2_pos #AU
        self.y1_ini = y1_pos #AU
        self.y2_ini = y2_pos #AU
        self.vx1_ini = vx1 #AU/yr
        self.vx2_ini = vx2 #AU/yr
        self.vy1_ini = vy1 #AU/yr
        self.vy2_ini = vy2 #AU/yr

    def solver_RK2(self,
                nsteps : int=100,
                delta_t : float=1./365.24):
        
        """
        It is a standard RK2 procedure.

        Parameters:
        Object1 : KeplerObject
            First object in the system
        Object2 : KeplerObject
            Second object in the system
        nsteps: int
            Number of time steps over which to solve the equations of motion

        Returns:
        Object1 : KeplerObject
            First object of the system with updated positions and velocities
        Object2 : KeplerObject
            Second object of the system with updated positions and velocities
        """
        vector = [
            self.Object1.x_pos, self.Object1.y_pos,
            self.Object2.x_pos, self.Object2.y_pos,
            self.Object1.vx, self.Object1.vy,
            self.Object2.vx, self.Object2.vy
        ]
        vector_new = np.array(vector)

        for i in range(nsteps):
            vector_inter = self.calculate_step(first_object=self.Object1,
                                               second_object=self.Object2,
                                               vector=vector_new,
                                               delta_t=delta_t/2.)
            
            vector_new = self.calculate_step(first_object=self.Object1,
                                             second_object=self.Object2,
                                             vector=vector_inter,
                                             delta_t=delta_t)
            
            self.Object1.update_vector(new_x=vector_new[0],
                                       new_y=vector_new[1],
                                       new_vx=vector_new[4],
                                       new_vy=vector_new[5])
            
            self.Object2.update_vector(new_x=vector_new[2],
                                       new_y=vector_new[3],
                                       new_vx=vector_new[6],
                                       new_vy=vector_new[7])

    def calculate_step(self,
                       first_object : KeplerObject,
                       second_object : KeplerObject,
                       vector : np.ndarray,
                       delta_t : float=1./365.24) -> np.ndarray:
        """
        Method allowing to update the position and velocities of the object as
        The method is a standard RK2 procedure.

        Parameters:
        Object1 : KeplerObject
            First object in the system
        Object2 : KeplerObject
            Second object in the system
        vector : np.ndarray
            Array containing the velocities and positions of both objects
        delta_t : float
            time step for time propagation

        Returns:
        vector_out : np.ndarray
            Updated vector field of the Kepler Equations after timestep delta_t
        """
        accx_1 = - first_object.acceleration(rel_pos=(first_object.x_pos - second_object.x_pos),
                                           mass_other=second_object.mass,
                                           x_pos_other=second_object.x_pos,
                                           y_pos_other=second_object.y_pos)
        accy_1 = - first_object.acceleration(rel_pos=(first_object.y_pos - second_object.y_pos),
                                           mass_other=second_object.mass,
                                           x_pos_other=second_object.x_pos,
                                           y_pos_other=second_object.y_pos)

        accx_2 = second_object.acceleration(rel_pos=(second_object.x_pos - first_object.x_pos),
                                              mass_other=first_object.mass,
                                              x_pos_other=first_object.x_pos,
                                              y_pos_other=first_object.y_pos)
        accy_2 = second_object.acceleration(rel_pos=(second_object.y_pos - first_object.y_pos),
                                              mass_other=first_object.mass,
                                              x_pos_other=first_object.x_pos,
                                              y_pos_other=first_object.y_pos)
        

        k = np.array([
            first_object.vx, first_object.vy,
            second_object.vx, second_object.vy,
            accx_1, accy_1,
            accx_2, accy_2
        ])

        vector_out = vector + delta_t * k
                
        return vector_out

    def solver_Verlet(self,
                nsteps : int=100,
                delta_t : float=1./365.24):
        
        f_object = self.Object1
        s_object = self.Object2
        vector_pos = np.array([f_object.x_pos, f_object.y_pos,
                               s_object.x_pos, s_object.y_pos])
        vector_v = np.array([f_object.vx, f_object.vy,
                             s_object.vx, s_object.vy])
        
        vector_pos_new = vector_pos
        vector_v_new = vector_v

        for i in range(nsteps):
            accx_1 = f_object.acceleration(rel_pos=(f_object.y_pos - s_object.y_pos),
                                            mass_other=s_object.mass,
                                            x_pos_other=s_object.x_pos,
                                            y_pos_other=s_object.y_pos)
            accy_1 = f_object.acceleration(rel_pos=(f_object.y_pos - s_object.y_pos),
                                            mass_other=s_object.mass,
                                            x_pos_other=s_object.x_pos,
                                            y_pos_other=s_object.y_pos)

            accx_2 = s_object.acceleration(rel_pos=(s_object.x_pos - f_object.x_pos),
                                                mass_other=f_object.mass,
                                                x_pos_other=f_object.x_pos,
                                                y_pos_other=f_object.y_pos)
            accy_2 = s_object.acceleration(rel_pos=(s_object.y_pos - f_object.y_pos),
                                                mass_other=f_object.mass,
                                                x_pos_other=f_object.x_pos,
                                                y_pos_other=f_object.y_pos)
            
            vector_a = np.array([accx_1, accy_1, accx_2, accy_2])

            vector_pos_new = vector_pos_new + vector_v_new * delta_t + 0.5 * vector_a * delta_t**2.
            rel_pos = vector_pos_new[0] - vector_pos_new[2]
            accx_1_new = - f_object.acceleration(rel_pos=rel_pos,
                                            mass_other=s_object.mass,
                                            x_pos_other=vector_pos_new[2],
                                            y_pos_other=vector_pos_new[3])
            rel_pos = vector_pos_new[1] - vector_pos_new[-1]
            accy_1_new = - f_object.acceleration(rel_pos=rel_pos,
                                            mass_other=s_object.mass,
                                            x_pos_other=vector_pos_new[2],
                                            y_pos_other=vector_pos_new[3])
            rel_pos = vector_pos_new[2] - vector_pos_new[0]
            accx_2_new = s_object.acceleration(rel_pos=rel_pos,
                                                mass_other=f_object.mass,
                                                x_pos_other=vector_pos_new[0],
                                                y_pos_other=vector_pos_new[1])
            rel_pos = vector_pos_new[-1] - vector_pos_new[1]
            accy_2_new = s_object.acceleration(rel_pos=rel_pos,
                                                mass_other=f_object.mass,
                                                x_pos_other=vector_pos_new[0],
                                                y_pos_other=vector_pos_new[1])

            vector_a_new = np.array([accx_1_new, accy_1_new, 
                                     accx_2_new, accy_2_new])
            
            vector_v_new = vector_v_new + 0.5 * (vector_a + vector_a_new) * delta_t

            vector_new = list(vector_pos_new) + list(vector_v_new)

            self.Object1.update_vector(new_x=vector_new[0],
                                       new_y=vector_new[1],
                                       new_vx=vector_new[4],
                                       new_vy=vector_new[5])
            
            self.Object2.update_vector(new_x=vector_new[2],
                                       new_y=vector_new[3],
                                       new_vx=vector_new[6],
                                       new_vy=vector_new[7])

        

    # calculer x(i+1) avec x_i, v_i et a_i
    # calculer v(i+1) à partir de a_i+1

    def total_energy(self) -> float:
        total_kinetic_energy = self.Object1.kinetic_energy() + self.Object2.kinetic_energy()
        potential_energy = - self.Object1.G * self.Object1.mass * self.Object2.mass
        relative_distance = np.sqrt((self.Object1.x_pos - self.Object2.x_pos)**2. + (self.Object1.y_pos - self.Object2.y_pos)**2.)
        potential_energy = potential_energy / relative_distance
        return total_kinetic_energy + potential_energy

    def reinitialize_system(self):
        self.Object1 = KeplerObject(mass=self.Object1.mass,
                                    x_pos=self.x1_ini,
                                    y_pos=self.y1_ini,
                                    vx=self.vx1_ini,
                                    vy=self.vy1_ini)
        self.Object2 = KeplerObject(mass=self.Object2.mass,
                                    x_pos=self.x2_ini,
                                    y_pos=self.y2_ini,
                                    vx=self.vx2_ini,
                                    vy=self.vy2_ini)
        

    
    
