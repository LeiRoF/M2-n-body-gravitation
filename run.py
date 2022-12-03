import os
from time import time

if not os.path.isdir("data"):
    os.mkdir("data")

print("\n⚙️ Compiling...")

if os.name == 'nt':
    os.system(r"powershell src\compile.ps1")
else:
    os.system("src/compile.sh")

print("✅ Done!")

print("\n🏃 Running...")
start = time()

if os.name == 'nt':
    os.system(".\data\main.exe")
else:
    os.system("./data/main.out")

end = time()
print("✅ Done!")

print(f"\n⌚ Execution time: {(end - start):.3} seconds")