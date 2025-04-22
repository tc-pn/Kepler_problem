"""
Visualization module.
"""
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
from SystemKeplerClass import KeplerSystem

def show_trajectory(System : KeplerSystem):
    traj_1 = np.array([list(coordinates[:2]) for coordinates in System.Object1.trajectory])
    traj_2 = np.array([list(coordinates[:2]) for coordinates in System.Object2.trajectory])

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

