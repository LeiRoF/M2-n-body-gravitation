import os
import numpy as np
import matplotlib.pyplot as plt
import regex as re
import multiprocessing

# Run a simulation if no data is found
if not os.path.isfile("data/speedup.txt"):
    import run

N = multiprocessing.cpu_count()

print("\n🔎 Retrieving results...")
speeds = np.empty(N)
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

print("\n⚙️ Generating plot...")
fig = plt.figure()
plt.subplot(121)
plt.plot(np.arange(N),speeds)
plt.title("Execution time")
plt.xlabel("Number of threads")
plt.ylabel("Time [s]")
plt.grid()

plt.subplot(122)
plt.plot(np.arange(N),speeds[0] / speeds)
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