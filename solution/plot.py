import os
from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation

N = 2000
steps = 1000
dt = 0.01

def load_file(filename, in_line): 
    data = np.fromfile(filename, dtype = np.float64) 
    print(data.shape)
    data = data.reshape(data.shape[0]//(in_line*3), in_line, 3) 
    return data 

energies = load_file("energy_fort_seq.dat", N)
bodies = load_file("position_fort_seq.bin", N)

print(energies.shape)
print(bodies.shape)

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
# 3D plot

ax = fig1.add_subplot(335, projection='3d')
title = ax.set_title('N-body gravitation')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

graph = ax.scatter(bodies[0,:,0], bodies[0,:,1], bodies[0,:,2], c = np.arange(N)/N, cmap = "hsv")

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