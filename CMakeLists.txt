cmake_minimum_required(VERSION 3.23)

project(NestorCFD C Fortran)

add_executable(nestor)

option(NANCHECK "Check for NaN values in certain routines" ON)

include(CheckSourceRuns)
set(CMAKE_REQUIRED_FLAGS "-march=native")
check_source_runs(C "
    #include <immintrin.h>
    int main() {
        __m512 a = _mm512_set1_ps(0); // AVX512F
        (void)a;
        return 0;
    }
" HOST_AVX512)

if(HOST_AVX512)
    set(AVX_DEFAULT ON)
else()
    set(AVX_DEFAULT OFF)
endif()

option(USE_AVX512 "Use AVX512 vector functions" ${HOST_AVX512})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(LDFLAGS_OPT -flto=auto -fwhole-program)
    set(FFLAGS_DBG -O0 -g -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing
	     -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic
		 -fbacktrace -fall-intrinsics)
    set(FFLAGS_OPT -g -O3 -march=native ${LDFLAGS_OPT})
    set(CFLAGS_DEF -O3 -g -Wall -Wextra -march=native ${LDFLAGS_OPT})
    set(OMP_OPT "-fopenmp")
else()
    set(LDFLAGS_OPT "-flto")
    set(FFLAGS_DBG "-O0 -g -fimplicit-none  -Wall -std=f2008")
    set(FFLAGS_OPT "-g -O3 -march=native ${LDFLAGS_OPT}")
    set(CFLAGS_DEF "-O3 -g -Wall -Wextra -march=native ${LDFLAGS_OPT}")
    set(OMP_OPT "-qopenmp")
endif()

if(USE_AVX512)
    list(APPEND FFLAGS_DBG "-D__USE_VINTRINSICS")
    list(APPEND FFLAGS_OPT "-D__USE_VINTRINSICS")
endif()

if (build_type STREQUAL debug)
  target_compile_options(nestor PRIVATE
    $<$<COMPILE_LANGUAGE:C>:${CFLAGS_DEF}>

    $<$<COMPILE_LANGUAGE:Fortran>:${FFLAGS_DBG}>    
  )
else()
  target_compile_options(nestor PRIVATE
    $<$<COMPILE_LANGUAGE:C>:${CFLAGS_DEF}>

    $<$<COMPILE_LANGUAGE:Fortran>:${FFLAGS_OPT}>    
  )
  target_link_options(nestor PRIVATE ${LDFLAGS_OPT})
endif()

set_property(SOURCE limiter.F90 res_sa.F90 viscosity.F90
             APPEND
             PROPERTY COMPILE_OPTIONS "${OMP_OPT}")

set(NESTOR_SOURCES )
list(APPEND NESTOR_SOURCES ad_inviscid_flux.f90
    ad_operators.f90
    ad_viscous_flux.F90
    algebraic_multigrid.f90
    bc_states.f90
    direct_solve.f90
    files.f90
    finterface.f90
    force.f90
    gauss_seidel.f90
    gcr.f90
    gradient.F90
    grid.f90
    grid_statistics.f90
    initialize.F90
    inputoutput.f90
    interface_jacobian.f90
    inviscid_flux.f90
    jacobian.f90
    limiter.F90
    linear_solver.f90
    line_implicit.f90
    lowlevel.f90
    lsq.F90
    messages.f90
    nestor.f90
    parameters.F90
    reorder.f90
    residual.F90
    res_sa.f90
    res_turb.f90
    ruge_stuben.f90
    sa_vars.f90
    solution.f90
    solution_vars.f90
    sort_routines.f90
    sparse_block_matrix.f90
    sparse_common.f90
    sparse_scalar_matrix.f90
    steady_solver.f90
    turb_bc.f90
    turb.f90
    utils.f90
    viscosity.f90
    viscous_flux.F90
    wall_distance.f90)

if(USE_AVX512)
    list(APPEND NESTOR_SOURCES vector_intrinsics.c vi_interface.f90)
endif()

target_sources(nestor
PRIVATE
    ${NESTOR_SOURCES}
    )
    
message(${NESTOR_SOURCES})
    
#    testing.f90
