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
CC = gcc
# Note: use "gfortran -O3" for best performance, but
#       don't use it until you're sure bugs are removed.
FFLAGS = -O0 -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing \
	     -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  	    \
		 -fbacktrace -fall-intrinsics -DNANCHECK
LDFLAGS= -flto=auto -fwhole-program
CFLAGS = -O3 -g -Wall -Wextra -march=native
# FFLAGS = -O2 -pg
#  FFLAGS = -g -pg -O3 -march=native $(LDFLAGS)
##########################################################
# VPATH = ..
##########################################################
# Check for vector intrinsics AVX512F and AVX512DQ
USE_AVX512 := $(shell grep -q "avx512f" /proc/cpuinfo && grep -q "avx512dq" /proc/cpuinfo && echo 1 || echo 0)
ifeq ($(USE_AVX512), 1)
FFLAGS += -D__USE_VINTRINSICS
endif
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
.SUFFIXES : .o .f90 .F90
.f90.o:
	$(FC) $(FFLAGS) -c $<

.c.o:
	$(CC) $(CFLAGS) -c $<

%.o: %.F90 # run c preprocessor
	$(FC) -cpp $(FFLAGS) -c $<

##########################################################
SDIR = .

OBCTS = $(SDIR)/lowlevel.o\
		$(SDIR)/messages.o\
		$(SDIR)/utils.o\
		$(SDIR)/parameters.o
ifeq ($(USE_AVX512), 1)
OBCTS +=$(SDIR)/vector_intrinsics.o\
		$(SDIR)/vi_interface.o
endif
OBCTS +=$(SDIR)/ad_operators.o\
		$(SDIR)/sort_routines.o\
		$(SDIR)/files.o\
		$(SDIR)/grid.o\
		$(SDIR)/reorder.o\
		$(SDIR)/grid_statistics.o\
		$(SDIR)/line_implicit.o\
		$(SDIR)/wall_distance.o\
		$(SDIR)/solution_vars.o\
		$(SDIR)/sa_vars.o\
		$(SDIR)/viscosity.o\
		$(SDIR)/turb.o\
		$(SDIR)/solution.o\
		$(SDIR)/sparse_common.o\
		$(SDIR)/inputoutput.o\
		$(SDIR)/initialize.o\
		$(SDIR)/turb_bc.o\
		$(SDIR)/sparse_block_matrix.o\
		$(SDIR)/sparse_scalar_matrix.o\
		$(SDIR)/direct_solve.o\
		$(SDIR)/bc_states.o\
		$(SDIR)/lsq.o\
		$(SDIR)/gradient.o\
		$(SDIR)/inviscid_flux.o\
		$(SDIR)/viscous_flux.o\
		$(SDIR)/ad_inviscid_flux.o\
		$(SDIR)/ad_viscous_flux.o\
		$(SDIR)/interface.o\
		$(SDIR)/limiter.o\
		$(SDIR)/res_sa.o\
		$(SDIR)/res_turb.o\
		$(SDIR)/residual.o\
		$(SDIR)/gauss_seidel.o\
		$(SDIR)/ruge_stuben.o\
		$(SDIR)/algebraic_multigrid.o\
		$(SDIR)/linear_solver.o\
		$(SDIR)/interface_jacobian.o\
		$(SDIR)/jacobian.o\
		$(SDIR)/force.o\
		$(SDIR)/gcr.o\
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

.PHONY:test
test: testing

testing: testing.o $(OBCTS)
	$(FC) $(FFLAGS) -o $@ $(SDIR)/testing.o $(filter-out $(SDIR)/nestor.o, $(OBCTS)) 
##########################################################
# Clean up
##########################################################
.PHONY:clean
clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.mod0
	rm -f nestor
