import os

if not os.path.isdir("data"):
        os.mkdir("data")

# TODO (add support for windows)
if os.name == 'nt':
    print("⚠️ Windows is not supported for this test")
    exit()

else:
    # Compile
    print("\n⚙️ Compiling...")
    os.system("src/compile.sh")
    print("✅ Done!")

    # Run bulk simulations
    print("\n🏃 Running simulation with different amount of threads...")
    os.system("./src/speedup.sh")
    print("✅ Done!")
