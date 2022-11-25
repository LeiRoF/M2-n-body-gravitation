import os
# os.system("gfortran -fno-backtrace .\main.f90 -o main.exe")
# os.system("main.exe")

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation

bodies = np.loadtxt("nbody.txt")
N, steps, dt = np.loadtxt("parameters.txt")
N = int(N)
steps = int(steps)

print(bodies.shape)

bodies = bodies.reshape(steps, N, 9)

print(bodies.shape)

print(bodies)

# fig = plt.figure()
# ax = fig.add_subplot(131,projection='3d')
# ax.set_title("Position of objects")
# ax.scatter(bodies[0, :,0], bodies[0, :,1], bodies[0, :,2])
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("z")
# ax.set_box_aspect((1,1,1))


# ax = fig.add_subplot(132,projection='3d')
# ax.set_title("Velocity on X repartion")
# ax.scatter(bodies[0, :,0], bodies[0, :,1], bodies[0, :,3])
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("vx")
# ax.set_box_aspect((1,1,1))

# ax = fig.add_subplot(133,projection='3d')
# ax.set_title("Velocity on Y repartion")
# ax.scatter(bodies[0, :,0], bodies[0, :,1], bodies[0, :,4])
# ax.set_xlabel("x")
# ax.set_ylabel("y")
# ax.set_zlabel("vy")
# ax.set_box_aspect((1,1,1))

# plt.show()


# ----------

def update_graph(t):
    graph._offsets3d = (bodies[t,:,0], bodies[t,:,1], bodies[t,:,2])
    title.set_text('Time={}'.format(t))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
title = ax.set_title('N-body gravitation')

graph = ax.scatter(bodies[0,:,0], bodies[0,:,1], bodies[0,:,2])

ani = matplotlib.animation.FuncAnimation(fig, update_graph, steps, 
                               interval=40, blit=False)

plt.show()