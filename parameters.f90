module common

    ! This module contains constant values that can be used by any of the other modules
    implicit none

    public
    !--------------------------------------------------------------------
    ! Constants
    integer , parameter :: p2 = selected_real_kind(P=15) !Double precision
    real(p2), parameter ::   zero = 0.0_p2
    real(p2), parameter ::   half = 0.5_p2
    real(p2), parameter ::    one = 1.0_p2
    real(p2), parameter ::    two = 2.0_p2
    real(p2), parameter ::  three = 3.0_p2
    real(p2), parameter ::   four = 4.0_p2
    real(p2), parameter ::   five = 5.0_p2
    real(p2), parameter ::    six = 6.0_p2
    real(p2), parameter ::  third = 1.0_p2/3.0_p2
    real(p2), parameter :: fourth = 1.0_p2/4.0_p2
    real(p2), parameter ::  sixth = 1.0_p2/6.0_p2
    real(p2), parameter :: my_eps = epsilon(one)  !Machine zero w.r.t. 1.0.
    real(p2), parameter :: pi = 3.141592653589793238_p2
    integer             :: ix = 1, iy = 2, iz = 3
    real(p2), parameter, dimension(5,5) :: canonical_array = reshape( (/1,0,0,0,0, 0,1,0,0,0, 0,0,1,0,0, 0,0,0,1,0, 0,0,0,0,1/), &
                                                            (/5,5/))    

    integer, parameter, dimension(3) :: version = (/0,1,1/) ! 0.1.1
end module common

module config
    ! this module loads/sets all config paramaters using a series of namelists.
    
    use common, only : p2, zero

    implicit none

    public

    !-------------------------------------------------------------------------
    ! PROJECT VARS (&project)
    character(80) ::      project_name = "default"    ! project name
    character(80) ::         grid_type = "ugrid"
    
    namelist / project / &
      project_name, grid_type


    !-------------------------------------------------------------------------
    ! INPUT/OUTPUT FILES (&inputoutput)
    logical       :: generate_tec_file_b = .true.  ! tecplot boundary file = T
    logical       :: generate_tec_file_v = .false. ! tecplot volume file   = F
    logical       :: write_data          = .false.
    logical       :: import_data         = .false.

    namelist / inputoutput / &
      generate_tec_file_b, generate_tec_file_v, &
      write_data         , import_data

    !-------------------------------------------------------------------------
    ! FREESTREAM CONDITIONS (&freestream)
    real(p2) ::    M_inf = 0.3_p2        ! Freestream Mach number
    real(p2) ::      aoa = 0.0_p2        ! Angle of attack in degrees (x --> z)
    real(p2) :: sideslip = 0.0_p2        ! Sideslip angle in degrees (x --> y)
    real(p2) :: gauge_pressure = 0.0_p2

    namelist / freestream / &
      M_inf, aoa, sideslip, gauge_pressure

    !-------------------------------------------------------------------------
    ! SOLVER SETTINGS (&solver)
    real(p2)               :: CFL                    = 0.5_p2
    logical                :: CFL_ramp               = .false.
    real(p2)               :: CFL_init               = 0.1_p2
    integer                :: CFL_start_iter         = 10
    integer                :: CFL_ramp_steps         = 100
    integer                :: solver_max_itr         = 1000
    real(p2)               :: solver_tolerance       = 1.0e-05_p2
    character(80)          :: method_inv_flux        = "roe"
    character(80)          :: method_inv_jac           = "roe"
    character(80)          :: solver_type            = "rk"
    character(80)          :: jacobian_method        = "analytical"
    real(p2), dimension(5) :: eig_limiting_factor    = (/ 0.1_p2, 0.1_p2, 0.1_p2, 0.1_p2, 0.1_p2 /)  !eigenvalue limiting factor
    real(p2), dimension(5) :: variable_ur            = (/ 1, 1, 1, 1, 1 /)  ! Variable under relaxation factors (only used in 
    logical                :: limit_update           = .false.
    logical                :: perturb_initial        = .false.
    logical                :: random_perturb         = .false.
    ! Closed loop method for limiting CFL in cells with large estimated change to prevent divergence
    
    namelist / solver / &
      CFL, CFL_ramp, CFL_init, CFL_start_iter, CFL_ramp, &
      solver_max_itr, solver_tolerance, &
      method_inv_flux, method_inv_jac, &
      solver_type, jacobian_method, eig_limiting_factor, &
      variable_ur, limit_update, perturb_initial

    !-------------------------------------------------------------------------
    ! AMG SETTINGS (&amg)
    logical                 :: use_amg              = .true.  
    character(80)           :: smoother             = "gs"    ! relaxation scheme type
    integer                 :: lrelax_sweeps        = 500     ! number of sweeps
    real(p2)                :: lrelax_tolerance     = 0.1_p2  ! relaxation tolerance (reduction)
    integer                 :: max_amg_levels       = 5
    
    namelist / amg / &
    use_amg, smoother, lrelax_sweeps, lrelax_tolerance, max_amg_levels

    !-------------------------------------------------------------------------
    ! GRADIENT SETTINGS (&gradient)
    character(80)           :: grad_method               = "lsq"
    integer                 :: accuracy_order       = 1
    character(80)           :: lsq_stencil          = "w_vertex"
    real(p2)                :: lsq_weight           = zero
    logical                 :: use_limiter          = .false.

    namelist / gradient / &
      grad_method, accuracy_order, lsq_stencil, use_limiter, lsq_weight

    contains
        
    subroutine read_nml_config(namelist_file)

        implicit none

        character(len=*), intent(in) :: namelist_file
        
        integer :: os

        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading the input file: ", namelist_file ,"...."
        write(*,*)

        open(unit=10,file=trim(namelist_file),form='formatted',status='old',iostat=os)
        
        ! check for file errors
        if (os /= 0) then ! open returned error
            write(*,*)
            write(*,*) "-------------------------------------------------------"
            write(*,*) "ERROR READING: ", namelist_file, ". CANNOT CONTINUE."
            write(*,*) "-------------------------------------------------------"
            stop
        endif

        ! Read config settings
        read(unit=10,nml=project)
        read(unit=10,nml=inputoutput)
        read(unit=10,nml=freestream)
        read(unit=10,nml=solver)
        read(unit=10,nml=amg)
        read(unit=10,nml=gradient)
        
    
        write(*,*)
        write(*,*) " List of given namelist variables and their values"
        write(*,*)

        write(*,*) "PROJECT VARIABLES (&project)"
        write(*,nml=project) 

        write(*,*)
        write(*,*) "INPUTOUTPUT SETTINGS (&inputoutput)"
        write(*,nml=inputoutput)
        
        write(*,*)
        write(*,*) "FREESTREAM SETTINGS (&freestream)"
        write(*,nml=freestream)

        write(*,*)
        write(*,*) "SOLVER SETTINGS (&solver)"
        write(*,nml=solver)
        
        write(*,*)
        write(*,*) "AMG SETTINGS (&amg)"
        write(*,nml=amg)

        write(*,*)
        write(*,*) "GRADIENT SETTINGS (&gradient)"
        write(*,nml=gradient)
        
        write(*,*)
        write(*,*) " End of Reading the input file: ",namelist_file,"..... "
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    end subroutine read_nml_config

end module config