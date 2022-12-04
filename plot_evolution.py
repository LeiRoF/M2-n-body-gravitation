import os
from copy import deepcopy
# os.system("gfortran -fno-backtrace -O3 .\main.f90 -o main.exe")
# os.system("main.exe")

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation

# Run a simulation if no data is found
if not os.path.isfile("data/parameters.txt"):
    import run

# ----------------------------------------------------------------------------------------------------
# Retrieving data

print("\n🔎 Retrieving results...")

N, steps, dt = np.loadtxt("data/parameters.txt")
N = int(N)
steps = int(steps)

positions = np.loadtxt("data/positions.txt")
velocities = np.loadtxt("data/velocities.txt")
accelerations = np.loadtxt("data/accelerations.txt")
energies = np.loadtxt("data/energies.txt")
barycenter = np.loadtxt("data/barycenter.txt")

positions = positions.reshape(steps, N, 3)
velocities = velocities.reshape(steps, N, 3)
accelerations = accelerations.reshape(steps, N, 3)

print("✅ Done!")

# ----------------------------------------------------------------------------------------------------
# Create plot

fig = plt.figure()

# ----------------------------------------------------------------------------------------------------
# All energies

ax_energy = fig.add_subplot(234)
ax_energy.set_title('Energies')
ax_energy.plot(np.arange(steps),energies[:,0],label="Potential")
ax_energy.legend()
ax_energy.plot(np.arange(steps),energies[:,1],label="Kinetic")
ax_energy.legend()
ax_energy.plot(np.arange(steps),energies[:,2],label="Total")
ax_energy.legend()

# ----------------------------------------------------------------------------------------------------
# 2D plots

ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)

ax2.set_title("2D projections")

ax1.set_aspect(1)
ax2.set_aspect(1)
ax3.set_aspect(1)

ax1.set_xlabel('x')
ax1.set_ylabel('z')

ax2.set_xlabel('x')
ax2.set_ylabel('y')

ax3.set_xlabel('y')
ax3.set_ylabel('z')

acc_xy = accelerations[:,:,0]**2 + accelerations[:,:,1]**2
acc_xz = accelerations[:,:,0]**2 + accelerations[:,:,2]**2
acc_yz = accelerations[:,:,1]**2 + accelerations[:,:,2]**2
graph1 = ax1.scatter(positions[0,:,0], positions[0,:,2], c = acc_xz[0,:], cmap = "gnuplot")
graph2 = ax2.scatter(positions[0,:,0], positions[0,:,1], c = acc_xy[0,:], cmap = "gnuplot")
graph3 = ax3.scatter(positions[0,:,1], positions[0,:,2], c = acc_yz[0,:], cmap = "gnuplot")

# ----------------------------------------------------------------------------------------------------
# Barycenter

ax_barycenter = fig.add_subplot(236, projection='3d')

ax_barycenter.set_xlabel('x')
ax_barycenter.set_ylabel('y')
ax_barycenter.set_zlabel('z')
ax_barycenter.set_title('Barycenter')
ax_barycenter.set_xlim3d([-1.0, 1.0])
ax_barycenter.set_ylim3d([-1.0, 1.0])
ax_barycenter.set_zlim3d([-1.0, 1.0])

graph_barycenter = ax_barycenter.scatter([barycenter[0, 0]], [barycenter[0, 1]], [barycenter[0, 2]])

# ----------------------------------------------------------------------------------------------------
# 3D plots

ax = fig.add_subplot(235, projection='3d')
title = ax.set_title('N-body gravitation')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

graph = ax.scatter(positions[0,:,0], positions[0,:,1], positions[0,:,2], c = np.sum(accelerations[0,:,:]**2,axis=-1), cmap = "gnuplot")

# ----------------------------------------------------------------------------------------------------
# Animation

def update_graph(t):
    graph._offsets3d = (positions[t,:,0], positions[t,:,1], positions[t,:,2])
    graph.set_array(np.sum(accelerations[t,:,:]**2,axis=-1))
    title.set_text('Time={}'.format(t))
    graph_barycenter._offsets3d = ([barycenter[t, 0]], [barycenter[t, 1]], [barycenter[t, 2]])
    
    tmp = deepcopy(positions[t,:,:2]) # x,y 
    graph2.set_offsets(tmp)
    graph2.set_array(acc_xy[t,:])

    tmp[:,1] = deepcopy(positions[t,:,2]) # x,z
    graph1.set_offsets(tmp)
    graph1.set_array(acc_xz[t,:])
    
    tmp[:,0] = deepcopy(positions[t,:,1]) # y,z
    graph3.set_offsets(tmp)
    graph3.set_array(acc_yz[t,:])

    return (graph, graph1, graph2, graph3), #, graph_barycenter

# ----------------------------------------------------------------------------------------------------
# Generate animation

print("\n🎬 Generating animation")
ani = matplotlib.animation.FuncAnimation(fig, update_graph, steps, interval=40, blit=False)
print("✅ Done")

# ----------------------------------------------------------------------------------------------------
# Plotting animation

print("\n🖼️ Showing it")
plt.show()
print("✅ Done")

# ----------------------------------------------------------------------------------------------------
# Saving animation

print("\n🎞️ Exporting it as mp4")
ani.save("data/animation.mp4")
print("✅ Done")