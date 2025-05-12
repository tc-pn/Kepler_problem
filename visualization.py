"""
Visualization module.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

from SystemKeplerClass import KeplerSystem

def show_trajectory_anim(System : KeplerSystem):
    traj_1 = np.array(System.Object1.positions)
    traj_2 = np.array(System.Object2.positions)

    fig, ax = plt.subplots(1,1)
    line, = ax.plot([], [], lw=2)
    #raise ValueError("")
    anim = FuncAnimation(fig, lambda x : animate(x, line=line, traj=traj_2), 
                         frames=len(traj_2),
                         interval=1000)
    #anim.save("test.gif")
    plt.show()
    pass

def animate(n : int, 
            line,
            traj : np.ndarray):
    line.set_xdata(traj[:n, 0])
    line.set_ydata(traj[:n, 1])
    return line,

def show_trajectory(System : KeplerSystem):
    traj_1 = np.array(System.Object1.positions)
    traj_2 = np.array(System.Object2.positions)

    fig, ax = plt.subplots(1,1, figsize=(8,8))
    ax.plot(traj_2[:,0], traj_2[:,1], 
            color="blue")
    ax.plot(traj_1[:,0], traj_1[:,1], 
            color="blue")
    ax.grid()
    plt.show(block=False)
    pass

def show_total_energy(System : KeplerSystem):
    
    pass
