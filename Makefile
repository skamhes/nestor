##########################################################
# Makefile for NESTOR
##########################################################
 PROGRAM = nestor
##########################################################
# Default path for Fortran Standard Library
# Installation instructions:
# https://github.com/fortran-lang/stdlib#getting-started
##########################################################
# LD_LIBRARY_PATH="/usr/local/include/fortran_stdlib/GNU-9.4.0/" -I$(LD_LIBRARY_PATH)
##########################################################
# MAKE VARIABLES
FC = gfortran
# Note: use "gfortran -O3" for best performance, but
#       don't use it until you're sure bugs are removed.
FFLAGS = -O0 -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing \
	     -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  	    \
		 -fbacktrace -fall-intrinsics
# FFLAGS = -O2 -pg -c 
# FFLAGS = -O3 -c 
##########################################################
# Suffix Rule for f90
# The first line says to make sure that each object file
# is newer than its corresponding .mod file
# The second line defines the list of suffixes 
# The third line says that to create the target (object 
# file) the corresponding command should be run on the 
# prerequisit (.f90 file)
# note: pattern rules are not considered for the default
# target which is why "all" is run when we run "make"
##########################################################
%.o: %.mod
.SUFFIXES : .o .f90
.f90.o:
	@echo "suffix"
	$(FC) $(FFLAGS) -c $<

##########################################################
SDIR = .

OBCTS = $(SDIR)/parameters.o\
		$(SDIR)/inputoutput.o\
		$(SDIR)/grid.o\
		$(SDIR)/solution.o\
		$(SDIR)/direct_solve.o\
		$(SDIR)/bc_states.o\
		$(SDIR)/lsq.o\
		$(SDIR)/gradient.o\
		$(SDIR)/inviscid_flux.o\
		$(SDIR)/interface.o\
		$(SDIR)/limiter.o\
		$(SDIR)/residual.o\
		$(SDIR)/steady_solver.o\
        $(SDIR)/nestor.o
##########################################################
# Make all
# This target doesnt actually get used it just exists to 
# follow common make conventions
all:$(PROGRAM)		
##########################################################
# Make executable "mg" 
# Note: use "gfortran -O2" for best performance, but
#       don't use it until you're sure bugs are removed.
##########################################################
$(PROGRAM): $(OBCTS)
	$(FC) $(FFLAGS) -o $@ $(OBCTS)

##########################################################
# Clean up
##########################################################
.PHONY:clean
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.mod0
	rm -f kcfd
