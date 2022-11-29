import os
from copy import deepcopy
os.system("gfortran -fno-backtrace -O3 .\main.f90 -o main.exe")
os.system("main.exe")

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation

energies = np.loadtxt("energies.txt")
bodies = np.loadtxt("nbody.txt")
barycenter = np.loadtxt("barycenter.txt")

N, steps, dt = np.loadtxt("parameters.txt")
N = int(N)
steps = int(steps)

bodies = bodies.reshape(steps, N, 9)

def update_graph1(t):
    tmp = deepcopy(bodies[t,:,:2]) # x,y 
    graph2.set_offsets(tmp)

    tmp[:,1] = deepcopy(bodies[t,:,2]) # x,z
    graph1.set_offsets(tmp)
    
    tmp[:,0] = deepcopy(bodies[t,:,1]) # y,z
    graph3.set_offsets(tmp)

    return graph1,


fig1 = plt.figure()

# ----------------------------------------------------------------------------------------------------
# Energies

ax_energy = fig1.add_subplot(333)
ax_energy.plot(np.arange(steps),energies[:,1],label="Potential")
ax_energy.legend()
ax_energy = fig1.add_subplot(336)
ax_energy.plot(np.arange(steps),energies[:,2],label="Kinetic")
ax_energy.legend()
ax_energy = fig1.add_subplot(339)
ax_energy.plot(np.arange(steps),energies[:,3],label="Total")
ax_energy.legend()

# ----------------------------------------------------------------------------------------------------
# All energies

ax_energy = fig1.add_subplot(338)
ax_energy.plot(np.arange(steps),energies[:,1],label="Potential")
ax_energy.legend()
ax_energy.plot(np.arange(steps),energies[:,2],label="Kinetic")
ax_energy.legend()
ax_energy.plot(np.arange(steps),energies[:,3],label="Total")
ax_energy.legend()

# ----------------------------------------------------------------------------------------------------
# 2D plots

ax1 = fig1.add_subplot(331)
ax2 = fig1.add_subplot(334)
ax3 = fig1.add_subplot(337)

ax1.set_aspect(1)
ax2.set_aspect(1)
ax3.set_aspect(1)

ax1.set_xlabel('x')
ax1.set_ylabel('z')

ax2.set_xlabel('x')
ax2.set_ylabel('y')

ax3.set_xlabel('y')
ax3.set_ylabel('z')

graph1 = ax1.scatter(bodies[0,:,0], bodies[0,:,2], c = np.arange(N)/N, cmap = "hsv")
graph2 = ax2.scatter(bodies[0,:,0], bodies[0,:,1], c = np.arange(N)/N, cmap = "hsv")
graph3 = ax3.scatter(bodies[0,:,1], bodies[0,:,2], c = np.arange(N)/N, cmap = "hsv")

ani1 = matplotlib.animation.FuncAnimation(fig1, update_graph1, steps, interval=40, blit=False)

# ----------------------------------------------------------------------------------------------------
# 3D plots

ax = fig1.add_subplot(335, projection='3d')
title = ax.set_title('N-body gravitation')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

graph = ax.scatter(bodies[0,:,0], bodies[0,:,1], bodies[0,:,2], c = np.arange(N)/N, cmap = "hsv")

# ----------------------------------------------------------------------------------------------------
# Barycenter

# fig1 = plt.figure()
# ax_barycenter = fig1.add_subplot(236, projection='3d')

# ax_barycenter.set_xlabel('x')
# ax_barycenter.set_ylabel('y')
# ax_barycenter.set_zlabel('z')

# graph_barycenter = ax_barycenter.scatter([barycenter[0, 1]], [barycenter[0, 2]], [barycenter[0, 3]])

# ----------------------------------------------------------------------------------------------------
# Animation

def update_graph(t):
    graph._offsets3d = (bodies[t,:,0], bodies[t,:,1], bodies[t,:,2])
    title.set_text('Time={}'.format(t))
    # graph_barycenter._offsets3d = ([barycenter[t, 1]], [barycenter[t, 2]], [barycenter[t, 3]])
    
    tmp = deepcopy(bodies[t,:,:2]) # x,y 
    graph2.set_offsets(tmp)

    tmp[:,1] = deepcopy(bodies[t,:,2]) # x,z
    graph1.set_offsets(tmp)
    
    tmp[:,0] = deepcopy(bodies[t,:,1]) # y,z
    graph3.set_offsets(tmp)

    return (graph, graph1, graph2, graph3), #, graph_barycenter

# ----------

ani = matplotlib.animation.FuncAnimation(fig1, update_graph, steps, interval=40, blit=False)

ani.save("animation.gif")
plt.show()