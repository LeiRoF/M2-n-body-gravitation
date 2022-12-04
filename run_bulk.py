import os

# TODO (add support for windows)
if os.name == 'nt':
    print("⚠️ Windows is not supported for this test")
    exit()

if not os.path.isdir("data"):
    os.mkdir("data")

print("\n⚙️ Compiling...")

if os.name == 'nt':
    os.system(r"powershell src\compile.ps1")
else:
    os.system("src/compile.sh")

print("✅ Done!")

print("\n🏃 Running simulation with different amount of threads...")

if os.name == 'nt':
    os.system(r"powershell .\src\speedup.ps1")
else:
    os.system("./src/speedup.sh")

print("✅ Done!")
