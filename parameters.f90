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
    real(p2), parameter :: three_half = 3.0_p2/2.0_p2
    real(p2), parameter ::    two = 2.0_p2
    real(p2), parameter ::  three = 3.0_p2
    real(p2), parameter ::   four = 4.0_p2
    real(p2), parameter ::   five = 5.0_p2
    real(p2), parameter ::    six = 6.0_p2
    real(p2), parameter ::  third = 1.0_p2/3.0_p2
    real(p2), parameter :: two_third = 2.0_p2/3.0_p2
    real(p2), parameter :: fourth = 1.0_p2/4.0_p2
    real(p2), parameter :: four_third = 4.0_p2/3.0_p2
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
    character(100)::   second_namelist = "empty"      ! allows you to load a second namelist other than nestor.nml
    character(100)::           io_path = "./"         ! path to any input and output files (grid, bc, data, etc.)
    
    namelist / project / &
      project_name, grid_type, second_namelist, io_path


    !-------------------------------------------------------------------------
    ! INPUT/OUTPUT FILES (&inputoutput)
    logical       :: generate_tec_file_b = .true.  ! tecplot boundary file = T
    logical       :: generate_tec_file_v = .false. ! tecplot volume file   = F
    logical       :: write_data          = .false.
    logical       :: import_data         = .false.
    logical       :: lift                = .false.
    logical       :: drag                = .false.
    real(p2)      :: area_reference      = 1.0_p2
    character(80) :: time_method         = "remaining"

    namelist / inputoutput / &
      generate_tec_file_b, generate_tec_file_v, &
      write_data         , import_data,         &
      lift, drag, area_reference, time_method

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
    real(p2)               :: CFL_max                = 1e+010
    real(p2)               :: CFL_min                = 1e-012
    integer                :: solver_max_itr         = 1000
    real(p2)               :: solver_tolerance       = 1.0e-05_p2
    character(80)          :: method_inv_flux        = "roe"
    character(80)          :: method_inv_jac         = "roe"
    character(80)          :: solver_type            = "rk"
    character(80)          :: jacobian_method        = "analytical"
    real(p2), dimension(5) :: eig_limiting_factor    = (/ 0.1_p2, 0.1_p2, 0.1_p2, 0.1_p2, 0.1_p2 /)  !eigenvalue limiting factor
    real(p2), dimension(5) :: variable_ur            = (/ 1, 1, 1, 1, 1 /)  ! Variable under relaxation factors (only used in 
    logical                :: limit_update           = .false.
    logical                :: perturb_initial        = .false.
    logical                :: random_perturb         = .false.
    logical                :: high_ar_correction     = .true.
    integer                :: gcr_max_projections    = 5
    real(p2)               :: gcr_reduction_target   = 0.9
    ! Closed loop method for limiting CFL in cells with large estimated change to prevent divergence
    
    namelist / solver / &
      CFL, CFL_ramp, CFL_init, CFL_start_iter, CFL_ramp, &
      solver_max_itr, solver_tolerance, &
      method_inv_flux, method_inv_jac, &
      solver_type, jacobian_method, eig_limiting_factor, &
      variable_ur, limit_update, perturb_initial, high_ar_correction, &
      gcr_max_projections, gcr_reduction_target

    !-------------------------------------------------------------------------
    ! AMG SETTINGS (&amg)
    logical                 :: use_amg              = .true.  
    character(80)           :: smoother             = "gs"    ! relaxation scheme type
    integer                 :: lrelax_sweeps        = 500     ! number of sweeps
    integer                 :: pre_sweeps           = 0       ! number of sweeps before AMG restriction
    integer                 :: post_sweeps          = 2       ! number of sweeps after AMG prolongation
    real(p2)                :: lrelax_tolerance     = 0.1_p2  ! relaxation tolerance (reduction)
    integer                 :: max_amg_levels       = 8
    integer                 :: min_amg_blcoks       = 1       ! minimum number of blocks before AMG will not further restrict
    character(1)            :: amg_cycle            = 'f'     ! amg cycle type.
    integer                 :: max_amg_cycles       = 16      ! Total complete AMG Cycles
    
    namelist / amg / &
    use_amg, smoother, lrelax_sweeps, lrelax_tolerance, max_amg_levels, pre_sweeps

    !-------------------------------------------------------------------------
    ! GRADIENT SETTINGS (&gradient)
    character(80)           :: grad_method               = "lsq"
    integer                 :: accuracy_order       = 1
    character(80)           :: lsq_stencil          = "w_vertex"
    real(p2)                :: lsq_weight           = zero
    logical                 :: use_limiter          = .false.

    namelist / gradient / &
      grad_method, accuracy_order, lsq_stencil, use_limiter, lsq_weight

    !-------------------------------------------------------------------------
    ! TURBULENCE SETTINGS (&turbulence)
      character(80)         :: turbulence_type       = 'inviscid'
      real(p2)              :: pr                    = 0.72_p2    ! Prandtl number for sea level air
      real(p2)              :: Re_inf                = 10000      ! Free stream reynolds number
      real(p2)              :: sutherland_constant   = 110.5_p2   ! (K) Sutherland's constant (C) for air
      real(p2)              :: ideal_gas_constant    = 287.058_p2 ! ideal gas constant for air (R)
      real(p2)              :: reference_temp        = 300        ! (K) T_inf in EQ 4.14.16 of I do Like CFD

    namelist / turbulence / &
      turbulence_type, pr, reference_temp, Re_inf, sutherland_constant, ideal_gas_constant


    !-------------------------------------------------------------------------
    ! DEBUG SETTINGS (&debug)
      integer :: gcr_verbosity = 0

    namelist / debug / &
      gcr_verbosity

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

        if (trim(second_namelist) /= 'empty' ) then ! I guess don't call your new namelist empty?
          ! Close the original namelist
          close(10)
          write(*,*) "Reading input from second namelist: ./", trim(second_namelist),"...."
          ! Open the new namelist with the same unit #
          open(unit=10,file=trim(second_namelist),form='formatted',status='old',iostat=os)
          if (os /= 0) then ! open returned error
            write(*,*)
            write(*,*) "-------------------------------------------------------"
            write(*,*) "ERROR READING: ", namelist_file, ". CANNOT CONTINUE."
            write(*,*) "-------------------------------------------------------"
            stop
          endif
          read(unit=10,nml=project)
        endif
        
        read(unit=10,nml=inputoutput,iostat=os)
        call nml_read_error_check(os,'I/O')
        rewind(10)

        read(unit=10,nml=freestream,iostat=os)
        call nml_read_error_check(os,'FREESTREAM')
        rewind(10)

        read(unit=10,nml=solver,iostat=os)
        call nml_read_error_check(os,'SOLVER')
        rewind(10)

        read(unit=10,nml=amg,iostat=os)
        call nml_read_error_check(os,'AMG')
        rewind(10)

        read(unit=10,nml=gradient,iostat=os)
        call nml_read_error_check(os,'GRADIENT')
        rewind(10)

        read(unit=10,nml=turbulence,iostat=os)
        call nml_read_error_check(os,'TURBULENCE')
        rewind(10)

        read(unit=10,nml=debug,iostat=os)
        call nml_read_error_check(os,'DEBUG')
        rewind(10)
        
    
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
        write(*,*) "TURBULENCE SETTINGS (&turbulence)"
        write(*,nml=turbulence)

        write(*,*)
        write(*,*) "VERBOSITY SETTINGS (&verbosity)"
        write(*,nml=debug)
        
        write(*,*)
        write(*,*) " End of Reading the input file: ",namelist_file,"..... "
        write(*,*) "-------------------------------------------------------"
        write(*,*)

        close(10)

        ! translate the char variables into integers
        call update_isettings

    end subroutine read_nml_config

    subroutine nml_read_error_check(status,nml_name)

      use iso_fortran_env, only : iostat_end, iostat_eor
      
      integer, intent(in)       :: status
      character(*), intent(in)  :: nml_name

      if (status == iostat_end .or. status == iostat_eor) then
        write(*,*) " NO ",nml_name," SETTINGS LOADED!!  USING DEFAULT SETTINGS."
      elseif(status /= 0) then
        write(*,*) " ERROR LOADING DEBUG SETTINGS!!"
      endif
    end subroutine nml_read_error_check

    subroutine update_isettings

      use utils

      implicit none

      call initialize_isettings

      ! update settings with imported namelist values

      select case(trim(time_method))
      case('remaining')
        itime_method = TM_REMAINING
      case('elapsed')
        itime_method = TM_ELAPSED
      case default
        write(*,*) ' time_method input "', trim(time_method),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select

      select case(trim(method_inv_flux))
      case('roe')
        imethod_inv_flux = IFLUX_ROE
      case default
        write(*,*) ' method_inv_flux input "', trim(method_inv_flux),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select
      
      select case(trim(method_inv_jac))
      case('roe')
        imethod_inv_jac = IJAC_ROE
      case default
        write(*,*) ' method_inv_jac input "', trim(method_inv_jac),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select
      
      select case(trim(solver_type))
      case('rk')
        isolver_type = SOLVER_RK
      case('explicit') 
        isolver_type = SOLVER_EXPLICIT
      case('implicit')
        isolver_type = SOLVER_IMPLICIT 
      case('gcr')
        isolver_type = SOLVER_GCR
      case default
        write(*,*) ' solver_type input "', trim(solver_type),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select
      
      select case(trim(jacobian_method))
      case('analytical')
        ijacobian_method = JAC_ANALYTIC
      case default
        write(*,*) ' jacobian_method input "', trim(jacobian_method),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select
      
      select case(trim(smoother))
      case('gs')
        ismoother = SMOOTH_GS
      case default
        write(*,*) ' smoother input "', trim(jacobian_method),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select

      select case(trim(grad_method))
      case('lsq')
        igrad_method = GRAD_LSQ
      case default
        write(*,*) ' grad_method input "', trim(jacobian_method),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select

      select case(trim(lsq_stencil))
      case('w_vertex')
        ilsq_stencil = LSQ_STENCIL_WVERTEX
      case default
        write(*,*) ' lsq_stencil input "', trim(jacobian_method),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select

      select case(trim(turbulence_type))
      case('inviscid')
        iturb_type = TURB_INVISCID
      case('laminar')
        iturb_type = TURB_LAMINAR
      case('rans')
        iturb_type = TURB_RANS
      case default
        write(*,*) ' turbulence_type input "', trim(jacobian_method),'" is invalid'
        write(*,*) ' error occured in update_isettings in utils.f90. Stopping...'
        stop
      end select
      
    end subroutine update_isettings

end module config