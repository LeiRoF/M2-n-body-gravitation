iterations = 10  # edit also the ./src/speedup.sh script

import os
import numpy as np
import matplotlib.pyplot as plt
import regex as re
import multiprocessing

N = multiprocessing.cpu_count() + 3

print("\n🔎 Retrieving results...")
times = np.loadtxt("data/speedup.txt")
print("✅ Done!")

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

print(times.shape)
print(times[0].shape)
print(threads)

print("\n⚙️ Generating plot...")
fig = plt.figure(figsize=(12,6))
plt.subplot(121)
for i in times:
    plt.plot(threads, i,"r", alpha=0.3)
plt.plot(threads, mean_times, "b")
plt.title("Execution time")
plt.xlabel("Number of threads")
plt.ylabel("Time [s]")
plt.grid()

plt.subplot(122)
for i in speedups:
    plt.plot(threads, i,"r", alpha=0.3)
plt.plot(threads, mean_speedups, "b")
plt.title("Speedup")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
plt.grid()
print("✅ Done!")

print("\n💾 Saving picture...")
plt.savefig("data/speedup.png")
print("✅ Done!")

print("\n📊 Plotting results...")
plt.show()
print("✅ Done!")