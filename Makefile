##########################################################
# Makefile for KCFD
##########################################################
 PROGRAM = nestor
##########################################################
# Default path for Fortran Standard Library
# Installation instructions:
# https://github.com/fortran-lang/stdlib#getting-started
##########################################################
# LD_LIBRARY_PATH="/usr/local/include/fortran_stdlib/GNU-9.4.0/" -I$(LD_LIBRARY_PATH)
##########################################################
# Suffix Rule for f90
# Note: use "gfortran -O2" for best performance, but
#       don't use it until you're sure bugs are removed.
##########################################################
.SUFFIXES : .o .f90
.f90.o:

	
	gfortran -O0 -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace -fall-intrinsics -c $<
#	gfortran -O2 -pg -c $<
#	gfortran -O3 -c $<
##########################################################
SDIR = .

OBCTS = $(SDIR)/parameters.o\
		$(SDIR)/inputoutput.o\
        $(SDIR)/nestor.o
##########################################################
# Make executable "mg" 
# Note: use "gfortran -O2" for best performance, but
#       don't use it until you're sure bugs are removed.
##########################################################
$(PROGRAM): $(OBCTS)
	gfortran -O0 -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace -fall-intrinsics -o $@ $(OBCTS)
#	gfortran -O2 -pg -o $@ $(OBCTS)
#	gfortran -O3 -o $@ $(OBCTS)
##########################################################
# Clean up
##########################################################
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.mod0
	rm -f kcfd
