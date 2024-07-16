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

        use config    , only : solver_type, accuracy_order, method_inv_flux, CFL, solver_max_itr, solver_tolerance, &
                            variable_ur, use_limiter, CFL_ramp, CFL_start_iter, CFL_ramp_steps, CFL_init
        
        use initialize, only : set_initial_solution

        use solution  , only : q, res, dtau, res_norm, res_norm_initial, lrelax_roc, lrelax_sweeps_actual, phi

        use grid      , only : cell, ncells

        use gradient  , only : init_gradients

        use residual  , only : compute_residual

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

        real(p2)                      :: CFL_multiplier, CFL_final, CFL_running_mult
        
        ! Formatting strings
        character(45) :: implicit_format = '(i10,6es13.3,a,i6,es12.1,i10.2,a,i2.2,es13.3)'
        character(35) :: explicit_format = '(i10,6es13.3,i10.2,a,i2.2,es13.3)'
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
        write(*,*) " inviscid_flux  = ", trim(method_inv_flux)
        

        if (CFL_ramp) then
            CFL_final = CFL       ! Final CFL after ramping
            CFL = CFL_init        ! Set the CFL to the initial CFL
            ! CFL multiplier to be performed after each step
            ! After each iter: CFL = CFL * CFL_mult
            CFL_multiplier = (CFL_final/CFL_init)**(one/CFL_ramp_steps)
            write(*,*) '    CFL Ramping = ENABLED'
            write(*,*) '    Initial CFL = ', CFL_init
            write(*,*) ' # of CFL Steps = ', CFL_ramp_steps
            write(*,*) ' CFL multiplier = ', CFL_multiplier
            write(*,*)
            write(*,*)
        else
            write(*,*) '    CFL Ramping = DISABLED'
            write(*,*) "            CFL = ", CFL
            write(*,*)
            write(*,*)
        endif

        if (accuracy_order == 2) then
            call init_gradients
        endif    

        ! Skipping importing data for now
        
        ! Print column headers
        if (trim(solver_type) == "implicit") then
            write(*,*) " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy       max-res", &
                "    |   proj     reduction       time    CFL"
            allocate( dq(5,ncells)) ! allocate du only if it needed
        else
            write(*,*) " Iteration   continuity   x-momemtum   y-momentum   z-momentum    energy       max-res    |   time     CFL"
        end if

        ! Initialize some miscellaneous variables
        lrelax_sweeps_actual = 0
        lrelax_roc = zero
        n_residual_evaluation = 0
        if (use_limiter) then
            allocate(phi(ncells))
        end if

        solver_loop : do while (i_iteration <= solver_max_itr)
            
            ! First compute the residual
            call compute_residual
            
            ! Compute residual norm
            call compute_residual_norm(res_norm)
            
            ! Iteration timer
            call dtime(values,time)
            
            ! Compute time remaining
            totalTime = time * real(solver_max_itr-i_iteration) ! total time remaining in seconds
            minutes = floor(totalTime/60.0)
            seconds = mod(int(totalTime),60)
            
            ! Allow the initial residual norm to increase for the first 5 iterations
            if ( i_iteration == 0 ) then
                res_norm_initial = res_norm
                minutes = 0
                seconds = 0
                do i = 1,5
                    ! Prevent res/res_norm_init = infinity
                    if (abs(res_norm_initial(i)) < 1e-016_p2) then
                        res_norm_initial(i) = one
                    end if
                end do
            elseif ( i_iteration <= 5 ) then
                do i = 1,5
                    if ( res_norm(i) > res_norm_initial(i) ) then
                        res_norm_initial(i) = res_norm(i)
                    end if
                end do
            endif

            ! Print out residual
            if ( trim(solver_type) == 'implicit' ) then
                write(*,implicit_format) i_iteration, res_norm(:), & 
                                         maxval(res_norm(:)/res_norm_initial(:)), &
                                         "   | ", lrelax_sweeps_actual, lrelax_roc, &
                                         minutes, ":", seconds, CFL
            else ! RK Explicit
                write(*,explicit_format) i_iteration, res_norm(:), &
                                         maxval(res_norm(:)/res_norm_initial(:)), &
                                         minutes, ":", seconds, CFL
            endif

            ! Check for convergence and exit if true
            if (maxval(res_norm(:)/res_norm_initial(:)) < solver_tolerance) then
                write(*,*) " Solution is converged!"
                exit solver_loop
            end if

            i_iteration = i_iteration + 1

            call compute_local_time_step_dtau

            ! March in pseudo-time to update u: u = u + du
            if (trim(solver_type) == "rk") then
                call explicit_pseudo_time_rk
            elseif (trim(solver_type) == 'explicit') then
                call explicit_pseudo_time_forward_euler
            elseif (trim(solver_type) == "implicit") then
                call implicit
            else
                write(*,*) " Unsopported iteration method: Solver = ", solver_type
            end if

            ! If using CFL ramp increase CFL
            if (CFL_ramp .and. (i_iteration < CFL_ramp_steps + CFL_start_iter) .and. i_iteration > CFL_start_iter) then
                CFL = CFL * CFL_multiplier
            elseif (CFL_ramp .and. (i_iteration == CFL_ramp_steps + CFL_start_iter)) then
                CFL = CFL_final
            end if

            ! check for stop file
            ! stop file can be created by typing "touch kcfdstop" in the working directory
            inquire (file = 'kcfdstop', exist = stop_me)
            if (stop_me) then
                write(*,*) "kcfdstop file found! Stopping iterations!"
                ! Delete the file that way it doesn't stop us next time.
                open(10,file = 'kcfdstop',status='old',iostat=ierr)
                if (ierr == 0) then
                    close(10,status ='delete',iostat = ierr) 
                    if (ierr == 0) then
                        write(*,*) 'kcfdstop successfully deleted!'
                    end if
                end if
                exit solver_loop
            end if

        end do solver_loop



    end subroutine steady_solve

    !********************************************************************************
    !
    ! This subroutine computes the residual L1 norm (average of the absolute value).
    !
    !********************************************************************************
    subroutine compute_residual_norm(res_norm)

        use common          , only : p2, zero
        use grid            , only : ncells
        use solution        , only : res, nq
    
        implicit none
    
        real(p2), dimension(nq), intent(out) :: res_norm
    
        !Local variables
        integer                :: i
    
        !Initialize the norm:
        res_norm(:) =  zero
    
        cell_loop : do i = 1, ncells
    
            res_norm(:) = res_norm(:) + abs( res(:,i) ) !L1 norm
    
        end do cell_loop
    
        res_norm(:) = res_norm(:) / real(ncells,p2)   !L1 norm
  
    end subroutine compute_residual_norm

    subroutine compute_local_time_step_dtau

        use common                  , only : half, p2
        use grid                    , only : ncells, cell
        use solution                , only : dtau, wsn
        use config                  , only : CFL

        implicit none

        integer :: i
        ! real(p2), dimension(ncells) :: viscous_dtau

        cell_loop : do i = 1,ncells
            dtau(i) = CFL * cell(i)%vol/( half * wsn(i) )

        end do cell_loop
    end subroutine compute_local_time_step_dtau

    subroutine explicit_pseudo_time_rk

    end subroutine explicit_pseudo_time_rk

    subroutine explicit_pseudo_time_forward_euler

        use common          , only : p2, half, one, zero

        use solution        , only : q, res, dtau, gammamo, gamma, gmoinv

        use grid            , only : cell, ncells

        use residual        , only : compute_residual

        use direct_solve    , only : gewp_solve

        real(p2), dimension(5) :: update_q
        integer i, os
        real(p2) :: H, rho_p, rho_T, theta, rho, uR2inv
        real(p2), dimension(5,5) :: preconditioner, pre_inv
        
        ! Compute the precondition matrix as described in https://doi.org/10.2514/3.12946
        ! Note: because we are not performing low-mach correction this is equivalent to the jacobian dW/dQ in equation (2)
        ! Eventually we will be adding low-mach preconditioning so this allows for future adaptability

        do i = 1,ncells
            ! test for low mach
            ! else
            uR2inv = one
            ! fi

            H = ((q(5,i))**2)*gmoinv + half * ( q(2,i)**2 + q(3,i)**2 + q(4,i)**2 )
            rho_p = gamma/q(5,i)
            rho_T = - (q(1,i)*gamma)/(q(5,i)**2)
            rho = q(1,i)*gamma/q(5,i)
            theta = (uR2inv) - rho_T*(gammamo)/(rho)
            
            preconditioner(1,:) = (/ theta,        zero,       zero,       zero,       rho_T                    /)
            preconditioner(2,:) = (/ theta*q(2,i), rho,        zero,       zero,       rho_T*q(2,i)             /)
            preconditioner(3,:) = (/ theta*q(3,i), zero,       rho,        zero,       rho_T*q(3,i)             /)
            preconditioner(4,:) = (/ theta*q(4,i), zero,       zero,       rho,        rho_T*q(4,i)             /)
            preconditioner(5,:) = (/ theta*H-one,  rho*q(2,i), rho*q(3,i), rho*q(4,i), rho_T*H + rho/(gamma-one)/)
            

            call gewp_solve(preconditioner, 5, pre_inv, os)
            if (os .ne. 0) then
                write(*,*) 'Error inverting precondition matrix at cell: ', i,' Stop!'
                stop
            end if
            
            update_q = - (dtau(i)/cell(i)%vol) * matmul(pre_inv,res(:,i))
            
            q(:,i) = q(:,i) + update_q
        end do
    end subroutine explicit_pseudo_time_forward_euler

    subroutine implicit

    end subroutine implicit

end module steady_solver