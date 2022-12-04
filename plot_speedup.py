iterations = 10  # edit also the ./src/speedup.sh script

import os
import numpy as np
import matplotlib.pyplot as plt
import regex as re
import multiprocessing

# Run a simulation if no data is found
if not os.path.isfile("data/speedup.txt"):
    import run

N = multiprocessing.cpu_count() + 3

print("\n🔎 Retrieving results...")
speeds = np.empty(N*iterations)
with open("data/speedup.txt") as f:

    i=0
    for line in f:
        res = re.search(r"real\s*([0-9])*m([0-9,.]*)s",line)
        if res is None:
            continue
        try:
            time = float(res.group(1)) * 60 + float(res.group(2))
        except: raise

        speeds[i] = time
        i+=1
        print(f"Execution took {time:3f} seconds on {i} threads")
print("✅ Done!")
speeds = speeds.reshape(N, iterations)

mean_speeds = np.empty(N)
std_speeds = np.empty(N)
for i in range(N):
    mean_speeds[i] = np.mean(speeds[:,i])
    std_speeds[i] = np.std(speeds[:,i])

speedups = np.empty_like(speeds)

for i in range(iterations):
    speedups[i,:] = mean_speeds[0]/speeds[i,:]

mean_speedups = np.empty(N)
std_speedups = np.empty(N)
for i in range(N):
    mean_speedups[i] = np.mean(speedups[:,i])
    std_speedups[i] = np.std(speedups[:,i])

threads = np.arange(1,N+1)

print("\n⚙️ Generating plot...")
fig = plt.figure()
# plt.subplot(121)
# for i in range(iterations):
#     plt.plot(threads,speeds[i,:],"r", alpha=0.3)
# plt.errorbar(threads, mean_speeds, std_speeds)
# plt.plot(threads, mean_speeds, "ob")
# plt.title("Execution time")
# plt.xlabel("Number of threads")
# plt.ylabel("Time [s]")
# plt.grid()

# plt.subplot(122)
for i in range(iterations):
    plt.plot(threads,speedups[i,:],"r", alpha=0.3)
# plt.errorbar(threads, mean_speedups, std_speedups)
# plt.plot(threads, mean_speedups, "ob")
plt.title("Speedup")
plt.xlabel("Number of threads")
plt.ylabel("Speedup")
# plt.grid()
print("✅ Done!")

print("\n💾 Saving picture...")
plt.savefig("data/speedup.png")
print("✅ Done!")

print("\n📊 Plotting results...")
plt.show()
print("✅ Done!")