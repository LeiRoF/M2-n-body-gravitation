THE CODE 
-------- 
Here is a sequential version of the N-body code. They are three files: 
    - main.f95: This is the main file that call subroutines from subroutine.f95 and settings from settings.f95 
    - subroutines.f95: Contains all the subroutines 
    - settings.f95: Contains the different settings of the program, including a subroutine that read the input.ini file 


INPUTS 
------ 
The file "input.ini", contains the value for the following parameters: 
    - npart: The number of particles 
    - nstep: The number of time steps 
    - dt: Time step 
    - epsilon: Gravitational softening 
    - init: Use (1) or not (0) the file "pos.dat" that contains pre-defined positions of 2000 particles. If used, the number of particles will be automatically set to 2000, whatever the value of npart. 

It allows to modify some parameters without the need to recompile the code. 


OUTPUTS 
------- 
They are two types of output files: 

    - "position_fort_seq.bin" contain the positions of each particle for each time step ("fort" is for "fortran" and "seq" is for "sequential").
        This is a binary file. To correctly read it in python, you can use the following function: 
  
        def load_file(filename, in_line): 
            data = np.fromfile(filename, dtype = np.float64) 
            data = data.reshape(data.shape[0]//(in_line*3), in_line, 3) 
            return data

        with filename="position_fort_seq.bin", and in_line the number of particles used for the simulation. 
  
    - "energy_fort_seq.dat" contains the potential and kinetic energy of the system for each time step. 1 line = 1 time step

The name of these two files can be modified in "settings.f95" 

  
COMPILATION 
----------- 
To compile the code, you can use the Makefile. Simply type "make" in your terminal. 
