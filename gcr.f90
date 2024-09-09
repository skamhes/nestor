module gcr

    use common , only : p2
    ! Method for computing the generalized conjugate residual, a Jacobian Free Newton-Krylov method.
    ! This is based on the paper: https://doi.org/10.2514/6.2021-0857 (A free version can be found on the NASA TRS).

    implicit none

    private 

    public :: gcr_run
    public :: gcr_failure_handler

    public :: GCR_SUCCESS,GCR_CFL_FREEZE, GCR_STALL, GCR_PRECOND_DIVERGE, GCR_PRECOND_STALL, GCR_REAL_FAIL
    
    integer, parameter :: GCR_SUCCESS = 0
    integer, parameter :: GCR_CFL_FREEZE = 1
    integer, parameter :: GCR_STALL   = 2
    integer, parameter :: GCR_PRECOND_DIVERGE = 3
    integer, parameter :: GCR_PRECOND_STALL = 4
    integer, parameter :: GCR_REAL_FAIL = 5

    real(p2) :: gcr_relaxation


    contains

    subroutine gcr_run(sol_update,iostat)

        use common      , only : p2

        use grid        , only : ncells

        use solution    , only : nq, jacobian_type, q

        implicit none

        real(p2),            dimension(nq,ncells), intent(out) :: sol_update
        integer ,                                  intent(out) :: iostat

        

        call gcr_solve(sol_update,iostat)

        if (iostat /= GCR_SUCCESS) return

        call gcr_real_check(q,sol_update,iostat)

        if (iostat /= GCR_SUCCESS) return

        call gcr_nl_control(sol_update,iostat)

    end subroutine gcr_run

    subroutine gcr_failure_handler(gcr_failure_code)

        use common , only : p2
        
        use config , only : CFL

        implicit none

        integer, intent(in) :: gcr_failure_code ! not actually used for now...

        CFL = CFL / 10.0_p2

    end subroutine gcr_failure_handler

    subroutine gcr_solve(gcr_final_update, iostat)

        use common      , only : p2, zero, one, half

        use config      , only : gcr_max_projections, gcr_reduction_target

        use grid        , only : ncells, cell

        use solution    , only : nq, inv_ncells, dtau, compute_primative_jacobian, jacobian_type, q, res, jac, &
                                 nl_reduction, n_projections

        use residual    , only : compute_residual

        use linear_solver, only: linear_relaxation, RELAX_FAIL_STALL, RELAX_FAIL_DIVERGE

        implicit none

        real(p2),            dimension(nq,ncells), intent(out) :: gcr_final_update   ! delta_Q_n+1
        integer,                                   intent(out) :: iostat                 ! status of gcr solve

        real(p2), dimension(nq,ncells,gcr_max_projections) :: delta_Q_k ! delta_Q_k
        real(p2), dimension(nq,ncells,gcr_max_projections) :: b_k   ! b_k

        real(p2)                                           :: b_k_mag ! | b_k |
        real(p2)                                           :: gcr_inner_prod           ! b_k^T * b_j = mu
        ! real(p2), dimension(nq,ncells)                     :: Q0    ! Initial primative vector at the start of outer iteration
        ! real(p2), dimension(nq,ncells)                     :: R0    ! Initial residual from the outer NL iteration
        real(p2)                                           :: RMS_R0
        real(p2)                                           :: update_rms ! || delta_Q_k ||
        real(p2)                                           :: update_mag ! l2 norm of preconditioned correction
        real(p2)                                           :: eps_frechet
        real(p2), dimension(nq,ncells)                     :: frechet_vector
        real(p2)                                           :: gamma_k ! gamma_k from EQ. 14
        real(p2), dimension(nq,ncells)                     :: gcr_residual
        real(p2)                                           :: gcr_res_rms
        
        real(p2), dimension(5,5)    :: preconditioner



        integer :: kdir, jdir ! projection direction indices
        integer :: icell
        integer :: ivar
        integer :: os

        ! Initialize the final correction
        gcr_final_update = zero
        gcr_residual = -res
        RMS_R0 = rms(nq,ncells,res,inv_ncells)
        
        project_loop : do kdir = 1,gcr_max_projections

            ! Run Preconditioner sweeps on the approximate jacobian
            call linear_relaxation(nq, jac, -gcr_residual, delta_Q_k(:,:,kdir),os)
            if (os == RELAX_FAIL_DIVERGE) then
                iostat = GCR_PRECOND_DIVERGE
                return
            elseif (os == RELAX_FAIL_STALL) then
                iostat = GCR_PRECOND_STALL
                return
            endif

            ! Compute Frechet Derivative
            update_rms = rms(nq,ncells,delta_Q_k(:,:,kdir),inv_ncells)
            update_mag = l2norm(nq,ncells,delta_Q_k(:,:,kdir))
            
            call compute_frechet( delta_Q_k(:,:,kdir),update_mag,update_rms,b_k(:,:,kdir) )

            ! Generate new search direction

            do icell = 1,ncells
                preconditioner = compute_primative_jacobian(q(:,icell))

                preconditioner = preconditioner * cell(icell)%vol/dtau(icell)

                b_k(:,icell,kdir) = b_k(:,icell,kdir) + matmul(preconditioner,delta_Q_k(:,icell,kdir))

            end do

            b_k_mag = l2norm(nq,ncells,b_k(:,:,kdir))

            ! Normalize the correction and search direction
            delta_Q_k(:,:,kdir) = delta_Q_k(:,:,kdir) / b_k_mag

            b_k(:,:,kdir)       = b_k(:,:,kdir)       / b_k_mag

            ! Normalize direction k to previous search direactions
            jdir = 1

            do while (jdir < kdir)
                ! Pretty sure this is just a modified Gram-Schmidt... suppose I could check...
                gcr_inner_prod = inner_product(ncells, b_k(:,:,kdir), b_k(:,:,jdir)) ! mu

                b_k(:,:,kdir) = b_k(:,:,kdir) - gcr_inner_prod * b_k(:,:,jdir)
                delta_Q_k(:,:,kdir)=delta_Q_k(:,:,kdir)-gcr_inner_prod * delta_Q_k(:,:,jdir)

                ! Renormalize
                b_k_mag = l2norm(nq,ncells,b_k(:,:,kdir)) 
                delta_Q_k(:,:,kdir) = delta_Q_k(:,:,kdir) / b_k_mag
                b_k(:,:,kdir)   = b_k(:,:,kdir)   / b_k_mag

            end do

            ! Compute projection of the current search direction and the previous residual
            gamma_k = inner_product(ncells, b_k(:,:,kdir), gcr_residual(:,:))

            ! if they are (nearly) orthogonal the solution has stalled and we should quit
            if ( gamma_k < 1.0e-03_p2 * l2norm(nq,ncells,gcr_residual) ) exit project_loop

            ! Update the final correction with the projection of delta_q_k
            gcr_final_update = gcr_final_update + gamma_k * delta_Q_k(:,:,kdir)

            ! Update the residual of the GCR system
            gcr_residual = gcr_residual - gamma_k * b_k(:,:,kdir)

            gcr_res_rms = rms(nq,ncells,gcr_residual,inv_ncells)

            nl_reduction = gcr_res_rms / RMS_R0

            if (nl_reduction < gcr_reduction_target) then
                iostat = GCR_SUCCESS
                n_projections = kdir
                return
            endif
        enddo project_loop

        ! We shouldn't be able to get here unless we ran out of projection directions
        iostat = GCR_STALL
        n_projections = gcr_max_projections

    end subroutine gcr_solve

    subroutine compute_frechet(sol_update,update_length,sol_rms,frechet_deriv)

        use common , only : p2, one

        use solution , only : nq, q, res

        use grid , only : ncells

        use residual , only : compute_residual

        implicit none

        real(p2), dimension(nq,ncells), intent(in) :: sol_update
        real(p2),                       intent(in) :: update_length
        real(p2),                       intent(in) :: sol_rms
        real(p2), dimension(nq,ncells), intent(out):: frechet_deriv

        real(p2), dimension(:,:), pointer :: q_n
        real(p2), dimension(:,:), pointer :: r_0
        real(p2)                          :: eps_frechet
        real(p2)                          :: frech_min_bound = 1.0e-07_p2

        ! compute eps to be used in the frechet derivative
        eps_frechet = max(sol_rms,one)*frech_min_bound

        ! move the lates solution vector to the temp vector q_n
        q_n => q
        r_0 => res

        ! Set q = q + eps * dq/|dq|
        ! I think nullifying and reallocating the solution and residual vectors should be faster than directly copying them to 
        ! temp variables.  At some point I may test that to confirm.
        nullify(q,res)

        allocate(  q(nq,ncells))
        allocate(res(nq,ncells))

        q = q_n + eps_frechet * sol_update / update_length

        call compute_residual

        frechet_deriv = update_length * ( res - r_0 ) / eps_frechet

        deallocate(  q)
        deallocate(res)

        q   => q_n
        res => r_0

        nullify(q_n,r_0)

    end subroutine compute_frechet

    subroutine gcr_real_check(sol_current,sol_update, iostat)

        ! Realizability chekc for the proposed nonlinear update from gcr_run

        use common  , only : p2,zero

        use solution, only : nq

        use grid    , only : ncells

        implicit none 

        real(p2), dimension(nq,ncells), intent(in) :: sol_current
        real(p2), dimension(nq,ncells), intent(in) :: sol_update
        integer                       , intent(out):: iostat

        integer :: icell
        
        do icell = 1,ncells
            if (sol_current(1,icell) + sol_update(1,icell) <= zero .OR. sol_current(5,icell) + sol_update(1,icell) <= zero) then
                iostat = GCR_REAL_FAIL
                return
            endif
        end do  

        iostat = GCR_SUCCESS

    end subroutine gcr_real_check

    subroutine gcr_nl_control(sol_update, iostat)

        ! Determine the optimum underrelaxation factor for the nonlinear solution update.

        use common , only : p2, zero, half, one, two

        use config , only : gcr_reduction_target

        use solution , only : nq, q, res, compute_primative_jacobian, dtau, inv_ncells

        use grid , only : ncells, cell

        use residual , only : compute_residual

        implicit none

        real(p2), dimension(nq,ncells), intent(in) :: sol_update
        integer,                        intent(out):: iostat

        real(p2), dimension(nq,ncells) :: q0, R0
        real(p2)                       :: residual_reduct_target
        real(p2)                       :: Rtau_rms, R0_rms
        real(p2)                       :: delQ_rms, Qn_rms
        real(p2)                       :: f_0, f_1, g_1     ! terms in the optimization equation 22
        real(p2), dimension(nq)        :: update_i, frechet_i, R0_i
        real(p2)                       :: eps_frechet
        real(p2)                       :: update_rms
        real(p2)                       :: ur_opt

        integer :: icell, ivar

        q0 = q
        r0 = res

        q = q0 + sol_update

        residual_reduct_target = half * (one + gcr_reduction_target)

        call compute_residual

        do icell = 1,ncells
            res(:,icell) = res(:,icell) +  &
                            matmul( compute_primative_jacobian(q0(:,icell)) * cell(icell)%vol/dtau(icell), sol_update(:,icell) )
        end do  

        Rtau_rms = rms(nq,ncells,res,inv_ncells)
        R0_rms   = rms(nq,ncells,R0 ,inv_ncells)

        if ( Rtau_rms / R0_rms < residual_reduct_target ) then
            ! q and res have already been updated and the residual has reduced. We can return successful
            iostat = GCR_SUCCESS
            return
        endif

        ! Check if the change is comperable to the computer percision
        delQ_rms = rms(nq,ncells,q,inv_ncells)
        Qn_rms   = rms(nq,ncells,q0,inv_ncells)

        if ( delQ_rms / Qn_rms < 1.0e-12_p2) then
            if ( Rtau_rms / R0_rms < one) then
                ! Any reduction can be considered a success at this point
                iostat = GCR_CFL_FREEZE
                return
            endif
        endif

        ! The residual did not reduce so now we will apply an underrelaxation factor to minimize the residual
        f_0 = R0_rms
        f_1 = Rtau_rms
        
        eps_frechet = max(one,Qn_rms) * 1.0e-07_p2
        update_rms = rms(nq,ncells,sol_update,inv_ncells)
        q = q0 + eps_frechet * sol_update / update_rms
        call compute_residual
        g_1 = zero

        do icell = 1,ncells
            update_i = matmul( compute_primative_jacobian(q0(:,icell)) * cell(icell)%vol/dtau(icell), sol_update(:,icell) )
            frechet_i = update_rms * (res(:,icell) - R0(:,icell)) / eps_frechet
            R0_i = R0(:,icell)
            do ivar = 1,nq
                g_1 = g_1 + (update_i(ivar) + frechet_i(ivar) + R0_i(ivar))**2
            end do
        end do
        g_1 = sqrt( g_1 * inv_ncells )

        ur_opt = ( g_1 - f_0 ) / (two * (f_1 - g_1) )

        q = q0 + ur_opt * sol_update

        ! Check convergence
        call compute_residual
        ! Using R0 since we're done with it and we want to leave the residual untouched
        do icell = 1,ncells
            R0(:,icell) = res(:,icell) + &
                        matmul( compute_primative_jacobian(q0(:,icell)) * cell(icell)%vol/dtau(icell), sol_update(:,icell) )
        end do  

        Rtau_rms = rms(nq,ncells,R0,inv_ncells)

        if ( Rtau_rms / R0_rms < residual_reduct_target .OR. delQ_rms / Qn_rms < 1.0e-12_p2) then
            ! q and res have already been updated and the residual has reduced.
            iostat = GCR_CFL_FREEZE
            return
        endif

        ! If we've made it this far we failed :(
        iostat = GCR_STALL

    end subroutine gcr_nl_control

    pure function rms(nq,ncells,vector,div)

        ! Function for computing the L2 norm of block vectors

        use common , only : p2,zero 

        implicit none

        integer, intent(in)                 :: nq, ncells
        real(p2),dimension(:,:), intent(in) :: vector
        real(p2),                intent(in) :: div     ! dividend
        real(p2)                            :: rms

        integer :: i,j

        rms = zero

        do i = 1,ncells
            do j = 1,nq
                rms = rms + vector(nq,ncells)**2
            end do
        end do

        rms = sqrt(rms*div)
    end function

    pure function l2norm(nq,ncells,vector)

        ! Function for computing the L2 norm of block vectors

        use common , only : p2,zero 

        implicit none

        integer, intent(in)                 :: nq, ncells
        real(p2),dimension(:,:), intent(in) :: vector
        real(p2)                            :: l2norm

        integer :: i,j

        l2norm = zero

        do i = 1,ncells
            do j = 1,nq
                l2norm = l2norm + vector(nq,ncells)**2
            end do
        end do

        l2norm = sqrt(l2norm)
    end function

    pure function inner_product(ncells,vector1,vector2)

        ! Function for computing the inner_product of block vectors

        use common , only : p2,zero 

        implicit none

        integer, intent(in)                 :: ncells
        real(p2),dimension(:,:), intent(in) :: vector1,vector2
        real(p2)                            :: inner_product

        integer :: i

        inner_product = zero

        do i = 1,ncells
            inner_product = inner_product + dot_product(vector1(:,i),vector2(:,i))
        end do

    end function

end module gcr