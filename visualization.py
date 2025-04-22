"""
Visualization module.
"""
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
from SystemKeplerClass import KeplerSystem

def show_trajectory(System : KeplerSystem):
    traj_1 = [coordinates[:1] for coordinates in System.Object1.trajectory]
    traj_2 = [coordinates[:1] for coordinates in System.Object2.trajectory]

    fig, ax = plt.subplots(1,1)
    line = ax.plot([], [], lw=2)
    anim = FuncAnimation(fig, lambda x : animate(x, line=line, traj=traj_2), frames=len(traj_2), interval=len(traj_2))
    plt.show(block=False)
    pass

def animate(n : int, 
            line,
            traj):
    line.set_xdata(traj[:n][0])
    line.set_ydata(traj[:n][1])
    return line

