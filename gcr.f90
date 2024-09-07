module gcr

    use common , only : p2
    ! Method for computing the generalized conjugate residual, a Jacobian Free Newton-Krylov method.
    ! This is based on the paper: https://doi.org/10.2514/6.2021-0857 (A free version can be found on the NASA TRS).

    implicit none

    private 

    public :: gcr_solve

    public :: GCR_SUCCESS, GCR_STALL, GCR_PRECOND_DIVERGE, GCR_PRECOND_STALL, GCR_REAL_FAIL
    
    integer, parameter :: GCR_SUCCESS = 0
    integer, parameter :: GCR_CFL_FREEZE = 1
    integer, parameter :: GCR_STALL   = 2
    integer, parameter :: GCR_PRECOND_DIVERGE = 3
    integer, parameter :: GCR_PRECOND_STALL = 4
    integer, parameter :: GCR_REAL_FAIL = 5

    real(p2) :: gcr_relaxation


    contains

    subroutine gcr_solve(gcr_final_correction, iostat)

        use common      , only : p2, zero, one, half

        use config      , only : gcr_max_projections, gcr_reduction_target

        use grid        , only : ncells, cell

        use solution    , only : res, nq, ndim, jac, Q, inv_ncells, gamma, gammamo, gmoinv, dtau, compute_primative_jacobian

        use jacobian    , only : compute_jacobian

        use linear_solver, only: linear_relaxation, RELAX_FAIL_STALL, RELAX_FAIL_DIVERGE

        real(p2), dimension(nq,ncells), intent(out)        :: gcr_final_correction   ! delta_Q_n+1
        integer,                        intent(out)        :: iostat                 ! status of gcr solve

        real(p2), dimension(nq,ncells,gcr_max_projections) :: gcr_precond_correction ! delta_Q_k
        real(p2), dimension(nq,ncells,gcr_max_projections) :: gcr_search_direction   ! b_k

        real(p2)                                           :: gcr_search_direction_mag ! | b_k |
        real(p2)                                           :: gcr_inner_prod           ! b_k^T * b_j = mu
        real(p2), dimension(nq,ncells)                     :: Q0    ! Initial primative vector at the start of outer iteration
        real(p2), dimension(nq,ncells)                     :: R0    ! Initial residual from the outer NL iteration
        real(p2)                                           :: RMS_R0
        real(p2)                                           :: correction_rms ! || delta_Q_k ||
        real(p2)                                           :: correction_mag ! l2 norm of preconditioned correction
        real(p2)                                           :: eps_frechet
        real(p2)                                           :: gcr_projection ! gamma_k from EQ. 14
        real(p2), dimension(nq,ncells)                     :: gcr_residual
        real(p2)                                           :: gcr_res_rms
        real(p2)                                           :: gcr_reduction
        
        real(p2), dimension(5,5)    :: preconditioner



        integer :: kdir, jdir ! projection direction indices
        integer :: icell
        integer :: ivar
        integer :: os

        ! Store the initial values.  For now this isn't terribly efficient.  Maybe I'll change it later...
        R0 = res
        RMS_R0 = rms(nq,ncells,R0,inv_ncells)
        gcr_residual = -res
        Q0 = Q ! This can probably be optimized a bit by juggling things around with pointer.  But for now we'll do it slightly
        ! less efficiently to ensure it's right.  Then we can optimize

        gcr_final_correction = zero


        ! Generate the initial Approximate Jacobian and compute an initial update (preconditioner)
        call compute_jacobian
        
        project_loop : do kdir = 1,gcr_max_projections
            call linear_relaxation(nq, jac, res, gcr_precond_correction(:,:,kdir),os)
            if (os == RELAX_FAIL_DIVERGE) then
                iostat = GCR_PRECOND_DIVERGE
                return
            elseif (os == RELAX_FAIL_STALL) then
                iostat = GCR_PRECOND_STALL
                return
            endif

            ! Compute Frechet Derivative
            correction_rms = zero
            frechet : do icell = 1,ncells
                do ivar = 1,nq
                    correction_rms = correction_rms + gcr_precond_correction(ivar,icell,kdir)**2
                end do
            end do frechet
            correction_rms = sqrt( correction_rms*inv_ncells )
            correction_mag = l2norm(nq,ncells,gcr_precond_correction(:,:,kdir))
            eps_frechet = max(one,correction_rms) * 1.0e-07_p2

            q(:,:) = q0(:,:) + eps_frechet * gcr_precond_correction(:,:,kdir) / correction_mag

            call compute_residual

            ! Generate new search direction

            gcr_search_direction_mag = zero 

            do icell = 1,ncells
                preconditioner = compute_primative_jacobian(q0(:,icell))

                preconditioner(:,:) = preconditioner(:,:) * cell(icell)%vol/dtau(icell)

                gcr_search_direction(:,icell,kdir) = matmul(preconditioner,gcr_precond_correction(:,icell,kdir)) + &
                                                     ( res(:,icell) - R0(:,icell) ) / eps_frechet
            end do

            gcr_search_direction_mag = l2norm(nq,ncells,gcr_search_direction(:,:,kdir))

            ! Normalize the correction and search direction
            gcr_precond_correction(:,:,kdir) = gcr_precond_correction(:,:,kdir) / gcr_search_direction_mag

            gcr_search_direction(:,:,kdir)   = gcr_search_direction(:,:,kdir)   / gcr_search_direction_mag

            ! Normalize direction k to previous search direactions
            jdir = 1

            do while (jdir < kdir)
                gcr_inner_prod = inner_product(ncells, gcr_search_direction(:,:,kdir), gcr_search_direction(:,:,jdir))

                gcr_search_direction(:,:,kdir) = gcr_search_direction(:,:,kdir) - gcr_inner_prod * gcr_search_direction(:,:,jdir)
                gcr_precond_correction(:,:,kdir)=gcr_precond_correction(:,:,kdir)-gcr_inner_prod * gcr_precond_correction(:,:,jdir)

                ! Renormalize
                gcr_search_direction_mag = l2norm(nq,ncells,gcr_search_direction(:,:,kdir)) 
                gcr_precond_correction(:,:,kdir) = gcr_precond_correction(:,:,kdir) / gcr_search_direction_mag
                gcr_search_direction(:,:,kdir)   = gcr_search_direction(:,:,kdir)   / gcr_search_direction_mag

            end do

            gcr_projection = inner_product(ncells, gcr_search_direction(:,:,kdir), gcr_residual(:,:))

            ! TODO: Add logic for the case where gcr_projection << mag(gcr_residual)
            if ( gcr_projection < 1.0e-03_p2 * l2norm(nq,ncells,gcr_residual) ) exit project_loop

            gcr_final_correction = gcr_final_correction + gcr_projection * gcr_precond_correction(:,:,kdir)

            gcr_residual = gcr_residual - gcr_projection * gcr_search_direction(:,:,kdir)

            gcr_res_rms = rms(nq,ncells,gcr_residual,inv_ncells)

            gcr_reduction = gcr_res_rms / RMS_R0

            if (gcr_reduction < gcr_reduction_target) then
                iostat = GCR_SUCCESS
                q = Q0
                res = R0
                return
            endif
        enddo project_loop

        iostat = GCR_STALL
        q = Q0
        res = R0

    end subroutine gcr_solve

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
            res = res + matmul( compute_primative_jacobian(q0(:,icell)) * cell(icell)%vol/dtau(icell), sol_update(:,icell) )
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
            R0 = res + matmul( compute_primative_jacobian(q0(:,icell)) * cell(icell)%vol/dtau(icell), sol_update(:,icell) )
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