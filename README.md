# M2-n-body-gravitation

## 🔌 Installation

First, you need to install by yourself:
- [Python3](https://www.python.org/downloads/)
- [Git](https://git-scm.com/downloads)
- [Gfortran](https://gcc.gnu.org/wiki/GFortranBinaries)

Then, you can clone the repository and install the dependencies:

Then:

1. Clone the repository

    ```bash
    git clone https://github.com/LeiRoF/M2-n-body-gravitation
    ```

2. Install the dependencies

    ```bash
    pip install -r requirements.txt
    ```

That's it!

## 🚀 Usage

0. Edit the `config.f90` file to put the parameters you want.

1. If you want to see the evolution (the result of the program):
   
   1. Run the program once, with all threads

        ```bash
        python run.py
        ```
    2. Plot the evolution
    
        ```bash
        python plot_evolution.py
        ```
2. If you want to see the speedup over threads:

    1. Run the program with all threads

        ```bash
        python run_bulk.py
        ```
    2. Plot the speedup
    
        ```bash
        python plot_speedup.py
        ```

