iterations = 10  # edit also the ./src/speedup.sh script

import os
import numpy as np
import matplotlib.pyplot as plt
import regex as re
import multiprocessing

N = multiprocessing.cpu_count() + 3

print("\nš Retrieving results...")
times = np.loadtxt("data/speedup.txt")
print("ā Done!")

mean_times = np.empty(N)

for i in range(times.shape[1]):
    mean_times[i] = np.mean(times[:,i])

speedups = np.copy(times)

for i in range(times.shape[0]):
    speedups[i] = times[i,0]/times[i]

mean_speedups = np.empty(N)
for i in range(times.shape[1]):
    mean_speedups[i] = np.mean(speedups[:,i])

threads = np.arange(1,N+1)

def amdahl(p, n):
    return 1/(1-p+(p/n))

print("\nāļø Generating plot...")
fig = plt.figure(figsize=(12,6))
plt.subplot(121)
for i in times:
    plt.plot(threads, i,"r", alpha=0.3)
plt.plot(threads, mean_times, "b", label="Mean time")
plt.legend()
plt.title("Execution time")
plt.xlabel("Number of threads")
plt.ylabel("Time [s]")
plt.grid()

plt.subplot(122)
for i in speedups:
    plt.plot(threads, i,"r", alpha=0.3)
plt.plot(threads, mean_speedups, "b", label="Mean speedup")
plt.plot(threads, amdahl(0.579, threads), "g--", label="Amdahl's law p=0.579")
plt.plot(threads, amdahl(0.231, threads), "y--", label="Amdahl's law p=0.231")
plt.legend()
plt.title("Speedup")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
plt.grid()
print("ā Done!")

print("\nš¾ Saving picture...")
plt.savefig("data/speedup.png")
print("ā Done!")

print("\nš Plotting results...")
plt.show()
print("ā Done!")