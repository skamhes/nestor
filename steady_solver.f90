module steady_solver

    use common , only : p2
    
    use grid   , only : ncells

    use solution , only : nq

    implicit none

    public :: steady_solve

    public :: compute_residual_norm
    public :: compute_local_time_step_dtau

    public :: dq ! solution update
    real(p2), dimension(:,:), pointer :: dq

    private
    
    real(p2), dimension(5,5) :: var_ur_array
    integer                  :: i_iteration ! I may make this public if I need it

    contains

    subroutine steady_solve
        use common    , only : p2, half, one, zero 

        ! use linear_solver , only :  lrelax_sweeps_actual, lrelax_roc

        use config    , only : solver_type, accuracy_order, inviscid_flux, CFL, solver_max_itr, solver_tolerance, &
                            variable_ur, use_limiter, CFL_ramp, CFL_start_iter, CFL_ramp_steps, CFL_init
        
        use initialize, only : set_initial_solution

        use solution  , only : q, res, dtau, res_norm, res_norm_initial

        use grid      , only : cell, ncells

        use gradient  , only : init_gradients

        ! use residual  , only : compute_residual

        implicit none

        integer                       :: i, n_residual_evaluation
        integer                       :: L1 = 1

        ! Timing Variables
        real                          :: time, totalTime
        real, dimension(2)            :: values
        integer                       :: minutes, seconds

        ! Stop file
        logical                       :: stop_me
        integer                       :: ierr

        real(p2)                      :: CFL_multiplier, CFL_end, CFL_running_mult
        
        ! Set explicit under-relaxation array
        var_ur_array = zero
        do i = 1,5
            var_ur_array(i,i) = variable_ur(i)
        end do

        ! Set initial solution (or import but we'll do that later...)
        call set_initial_solution

        i_iteration = 0

        write(*,*) " ---------------------------------------"
        write(*,*) " Begin Nonlinear iterations"
        write(*,*) 
        write(*,*)
        write(*,*) "    solver_type = ", trim(solver_type)
        write(*,'(a,i1)') " accuracy_order = ", accuracy_order
        write(*,*) " inviscid_flux  = ", trim(inviscid_flux)
        write(*,*) "            CFL = ", CFL

        call init_gradients
        



    end subroutine steady_solve

    subroutine compute_residual_norm

    end subroutine compute_residual_norm

    subroutine compute_local_time_step_dtau

    end subroutine compute_local_time_step_dtau

end module steady_solver