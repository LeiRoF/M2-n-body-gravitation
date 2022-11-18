import os
os.system("gfortran -fno-backtrace .\main.f90 -o main.exe")
os.system("main.exe")

import matplotlib.pyplot as plt
import numpy as np

bodies = np.loadtxt("nbody.txt")

print(bodies.shape)
print(bodies)

fig = plt.figure()
ax = fig.add_subplot(131,projection='3d')
ax.set_title("Position of objects")
ax.scatter(bodies[:,0], bodies[:,1], bodies[:,2])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_box_aspect((1,1,1))


ax = fig.add_subplot(132,projection='3d')
ax.set_title("Velocity on X repartion")
ax.scatter(bodies[:,0], bodies[:,1], bodies[:,3])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("vx")
ax.set_box_aspect((1,1,1))

ax = fig.add_subplot(133,projection='3d')
ax.set_title("Velocity on Y repartion")
ax.scatter(bodies[:,0], bodies[:,1], bodies[:,4])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("vy")
ax.set_box_aspect((1,1,1))

plt.show()