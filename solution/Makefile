# Choose compiler
FC = gfortran

# Options and Path
FFLAGS = -O3


##############################################################################
#----------------------------------------------------------------------------

#.MAKEOPTS: -k -s

.SUFFIXES: .f95

.f95.o:
	$(FC) $(FFLAGS) -c $*.f95



OBJ    = settings.o \
		 subroutines.o \

# List all the "rules" that must be executed when "make" is typed in the command line
ALL: Test PartialClean
	echo "!!!        Compilation OK        !!!"

# Rule to compile the main and all the modules
Test : $(OBJ) main.o
	$(FC) $(FFLAGS) main.o $(OBJ) -o main

# Rule to remove all the intermediate files
PartialClean:
	rm *.o *.mod

# Rule to remove all the intermediate files and the executable
clean:
	rm *.o BGM *.mod

##############################################################################
#----------------------------------------------------------------------------
