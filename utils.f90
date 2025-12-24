module utils

    ! This module is used to convert strings to integer constants and vice versa.
    ! This should speed up the time it takes for conditional loops to execute.
    ! All variables will be the same as their string equivalent with i (integer) prefixed to the beginning of the variable name.
    ! Defaults are 0, with the exception of boundary conditions which borrows from FUN3D for its values

    implicit none

    private

    public itime_method, TM_REMAINING, TM_ELAPSED
    public imethod_inv_flux, imethod_inv_jac, ientropy_fix, isolver_type, ijacobian_method
    public IFLUX_ROE, IFLUX_ROE_LM, IJAC_ROE, IJAC_ROE_LM, SOLVER_RK, SOLVER_EXPLICIT, SOLVER_IMPLICIT, SOLVER_GCR, JAC_ANALYTIC
    public IJAC_RHLL, IJAC_HLL, IJAC_RUSANOV, IHARTEN, IMAVRIPLIS
    public ismoother, SMOOTH_GS
    public igrad_method, ilsq_stencil, GRAD_LSQ, LSQ_STENCIL_WVERTEX, LSQ_STENCIL_NN
    public iturb_type, TURB_DES, TURB_INVISCID, TURB_LAMINAR, TURB_LES, TURB_RANS
    public ibc_type, BC_BACK_PRESSURE, BC_FARFIELD, BC_TANGENT, BC_VISC_STRONG

    public initialize_isettings

    !I/O
    integer :: itime_method

    integer, parameter :: TM_REMAINING = 0
    integer, parameter :: TM_ELAPSED   = 1

    ! Solver settings
    integer :: imethod_inv_flux, imethod_inv_jac, ientropy_fix, isolver_type, ijacobian_method

    integer, parameter :: IFLUX_ROE    = 0
    integer, parameter :: IFLUX_ROE_LM = 1
    
    integer, parameter :: IJAC_ROE     = 0
    integer, parameter :: IJAC_ROE_LM  = 1
    integer, parameter :: IJAC_RUSANOV = 2
    integer, parameter :: IJAC_HLL     = 3
    integer, parameter :: IJAC_RHLL    = 4

    integer, parameter :: IHARTEN      = 0
    integer, parameter :: IMAVRIPLIS   = 1
    
    integer, parameter :: SOLVER_RK       = 0
    integer, parameter :: SOLVER_EXPLICIT = 1
    integer, parameter :: SOLVER_IMPLICIT = 2
    integer, parameter :: SOLVER_GCR      = 3

    integer, parameter :: JAC_ANALYTIC = 0

    ! AMG settings
    integer :: ismoother

    integer, parameter :: SMOOTH_GS = 0 ! Gauss-Seidel

    ! integer :: iamg_cycle

    ! Gradient settings
    integer :: igrad_method
    integer :: ilsq_stencil

    integer, parameter :: GRAD_LSQ = 0

    integer, parameter :: LSQ_STENCIL_WVERTEX = 0
    integer, parameter :: LSQ_STENCIL_NN      = 1

    ! Turbulence settings
    integer :: iturb_type

    integer, parameter :: TURB_INVISCID = 0
    integer, parameter :: TURB_LAMINAR  = 1
    integer, parameter :: TURB_RANS     = 2
    integer, parameter :: TURB_DES      = 3
    integer, parameter :: TURB_LES      = 4

    ! Boundary Conditions
    integer, dimension(:), allocatable :: ibc_type !boundary conditions stored as integers (using the same numbering as FUN3D)

    integer, parameter :: BC_VISC_STRONG   = 4000
    integer, parameter :: BC_FARFIELD      = 5000
    integer, parameter :: BC_TANGENT       = 3000 ! used for symmetry planes and slip walls
    integer, parameter :: BC_BACK_PRESSURE = 5051

    contains

    subroutine initialize_isettings
        
        implicit none

        itime_method     = 0
        imethod_inv_flux = 0
        imethod_inv_jac  = 0
        isolver_type     = 0
        ijacobian_method = 0
        ismoother        = 0
        igrad_method     = 0
        ilsq_stencil     = 0
        iturb_type       = 0

    end subroutine initialize_isettings


end module utils