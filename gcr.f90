module gcr

    use common , only : p2

    use config , only : gcr_verbosity
    ! Method for computing the generalized conjugate residual, a Jacobian Free Newton-Krylov method.
    ! This is based on the paper: https://doi.org/10.2514/6.2021-0857 (A free version can be found on the NASA TRS).

    implicit none

    private 

    public :: gcr_run
    public :: gcr_CFL_control

    public :: GCR_SUCCESS,GCR_CFL_FREEZE, GCR_STALL, GCR_PRECOND_DIVERGE, GCR_PRECOND_STALL, GCR_REAL_FAIL
    
    integer, parameter :: GCR_SUCCESS = 0
    integer, parameter :: GCR_CFL_FREEZE = 1
    integer, parameter :: GCR_STALL   = 2
    integer, parameter :: GCR_PRECOND_DIVERGE = 3
    integer, parameter :: GCR_PRECOND_STALL = 4
    integer, parameter :: GCR_REAL_FAIL = 5
    integer, parameter :: GCR_PREAL_FAIL = 6

    interface clear_jacobian_arrays
        module procedure clear_jacobian_arrays_B
        module procedure clear_jacobian_arrays_S
    end interface clear_jacobian_arrays

    contains

    subroutine gcr_run(iostat)

        use common      , only : p2

        use grid        , only : ncells

        use solution_vars, only : nq, jacobian_type, q, res

        use turb        , only : nturb, turb_var, turb_res

        implicit none
        
        integer ,          intent(out) :: iostat

        real(p2), dimension(nq,ncells)    :: sol_update_f
        real(p2), dimension(ncells,nturb) :: sol_update_t
        real(p2)                          :: gcr_res_rms

        real(p2), dimension(:,:), pointer :: q_0, r_n ! Initial solution (q0) and res (r_n)
        real(p2), dimension(:,:), pointer :: t_0, tr_n ! Ditto for turb vars


        q_0 => q
        r_n => res

        if (associated(turb_var)) t_0  => turb_var
        if (associated(turb_res)) tr_n => turb_res
        
        call gcr_solve_scratch(sol_update_f,sol_update_t,gcr_res_rms,iostat)

        if (iostat /= GCR_SUCCESS) return

        call gcr_real_check(q,sol_update_f,iostat) ! for now we don't need to do a real check on the turb

        if (iostat /= GCR_SUCCESS) return

        call gcr_nl_control(sol_update_f,sol_update_t,gcr_res_rms,iostat)

    end subroutine gcr_run

    subroutine gcr_solve_scratch(gcr_final_update_f,gcr_final_update_t, gcr_res_rms, iostat)

        use common      , only : p2, zero, one

        use config      , only : gcr_max_projections, gcr_reduction_target, amg_cycle

        use grid        , only : ncells, cell

        use solution_vars    , only : nq, inv_ncells, jacobian_type, res, jac, &
                                 nl_reduction, n_projections, q
        
        use solution    , only : compute_primative_jacobian

        use residual    , only : compute_residual

        use linear_solver, only: linear_relaxation, RELAX_FAIL_STALL, RELAX_FAIL_DIVERGE, build_A_BCSM, multilevel_cycle, &
                                 build_Dinv_array

        use algebraic_multigird, only : convert_amg_c_to_i

        use turb , only : nturb, turb_jac, turb_res, turb_var

        use utils , only : iflow_type, FLOW_LAMINAR

        use lowlevel , only : merge_array

        implicit none

        real(p2), dimension(nq,ncells),    intent(out) :: gcr_final_update_f     ! delta_Q_n+1
        real(p2), dimension(ncells,nturb), intent(out) :: gcr_final_update_t     ! delta_Q_n+1
        real(p2),                          intent(out) :: gcr_res_rms            ! rms of the residual of converged update
        integer,                           intent(out) :: iostat                 ! status of gcr solve

        real(p2), dimension(nq,ncells)                        :: r_k                    ! gcr residual
        real(p2), dimension(ncells,nturb)                     :: r_k_t                  ! gcr residual
        real(p2), dimension(nq,ncells,gcr_max_projections)    :: b_k
        real(p2), dimension(ncells,nturb,gcr_max_projections) :: b_k_t
        real(p2), dimension(nq,ncells,gcr_max_projections)    :: dQ_k
        real(p2), dimension(ncells,nturb,gcr_max_projections) :: dQ_k_t
        real(p2)                                              :: norm_b_k_inv
        real(p2)                                              :: norm_r_k
        real(p2)                                              :: norm_dQ_k
        real(p2)                                              :: rms_r_k
        real(p2)                                              :: rms_r_0
        real(p2)                                              :: rms_Q_n
        real(p2)                                              :: mu, gamma_k ! inner products


        ! Variables for preconditioning matrix M

        real(p2), dimension(:,:,:), pointer     :: V
        real(p2), dimension(:,:),   pointer     :: Vt   ! Values (5x5 block matrix) plus corresponding index
        integer , dimension(:),     pointer     :: C, Ct   ! Column index of each value
        integer , dimension(:),     pointer     :: R, Rt   ! Start index of each new row
        integer                                 :: nnz, nnzt
        real(p2), dimension(:,:,:), pointer     :: Dinv
        real(p2), dimension(:,:),   pointer     :: Dinvt

        integer :: cycle_type

        integer :: iturb, kdir, jdir ! projection direction indices
        integer :: os

        ! Nullify pointers to avoid undefined behavior
        nullify(V,C,R,Dinv)
        nullify(Vt,Ct,Rt,Dinvt)

        ! Initialize some variables
        gcr_final_update_f = zero
        ! call merge_array(-res,-turb_res, ncells, nq, nturb, r_k)
        r_k              = -res
 
        ! Build M (A Approx) for precondition solve of flow variables
        allocate(R(ncells+1))
        allocate(Dinv(5,5,ncells))
        call build_A_BCSM(ncells,cell,jac,V,C,R,nnz=nnz)
        call build_Dinv_array(ncells,jac,Dinv)
        ! Build M for precondition solve of turbulent variables
        if (iflow_type > FLOW_LAMINAR) then
            allocate(Rt(ncells+1))
            allocate(Dinvt(ncells,nturb))
            call build_A_BCSM(ncells,cell,turb_jac,Vt,Ct,Rt,nturb,nnz=nnzt)
            call build_Dinv_array(ncells,turb_jac,nturb,Dinvt)
            r_k_t = -turb_res
            gcr_final_update_t = zero
        end if
        
        rms_r_0          = rms(nq,nturb,ncells,r_k,r_k_t   ,inv_ncells)
        rms_Q_n          = rms(nq,nturb,ncells,q  ,turb_var,inv_ncells)


        ! Solve the preconditioner
        cycle_type = convert_amg_c_to_i(amg_cycle)

        proj_loop : do kdir = 1,gcr_max_projections
            
            ! Compute correction for flow
            ! keep_A = .true. so that V,C,R, and Dinv do not have to be rebuilt
            call multilevel_cycle(ncells,nq, V, C, R, -r_k, Dinv,cycle_type,.true.,dQ_k(:,:,kdir),os)
            ! Compute correction for turbulence eqs
            do iturb = 1,nturb ! nturb is set to zero for laminar/inviscid flow
                call multilevel_cycle(ncells, Vt(:,iturb), Ct, Rt, -r_k_t(:,iturb), Dinvt(:,iturb), cycle_type, &
                                      .true., dQ_k_t(:,iturb,kdir), os)
            end do


            if (os == RELAX_FAIL_DIVERGE) then
                iostat = GCR_PRECOND_DIVERGE
                ! Clear Jacobian arrays
                call clear_jacobian_arrays(V,C,R,Dinv)
                call clear_jacobian_arrays(Vt,Ct,Rt,Dinvt)
                return
            elseif (os == RELAX_FAIL_STALL) then
                iostat = GCR_PRECOND_STALL
                ! Clear Jacobian arrays
                call clear_jacobian_arrays(V,C,R,Dinv)
                call clear_jacobian_arrays(Vt,Ct,Rt,Dinvt)
                return
            endif
            
            ! Generate the new search direction
            norm_dQ_k = l2norm(nq,nturb,ncells,dQ_k(:,:,kdir),dQ_k_t(:,:,kdir))

            call compute_frechet(q,res,turb_var,turb_res, dQ_k(:,:,kdir),dQ_k_t(:,:,kdir),norm_dQ_k,rms_Q_n, &
                                    b_k(:,:,kdir),b_k_t(:,:,kdir),os)

            if (os == GCR_PREAL_FAIL) then
                iostat = GCR_PREAL_FAIL
                ! Clear Jacobian arrays
                call clear_jacobian_arrays(V,C,R,Dinv)
                call clear_jacobian_arrays(Vt,Ct,Rt,Dinvt)
                return
            endif

            ! Orthonormalize
            norm_b_k_inv = one / l2norm(nq,nturb,ncells,b_k(:,:,kdir),b_k_t(:,:,kdir))
            
            b_k(  :,:,kdir) = b_k(  :,:,kdir) * norm_b_k_inv
            dQ_k(  :,:,kdir) = dQ_k(  :,:,kdir) * norm_b_k_inv
            if (iflow_type > FLOW_LAMINAR) then
                b_k_t(:,:,kdir) = b_k_t(:,:,kdir) * norm_b_k_inv
                dQ_k_t(:,:,kdir) = dQ_k_t(:,:,kdir) * norm_b_k_inv
            endif
            do jdir = 1,kdir - 1
                mu = inner_product(nq,nturb,ncells,b_k(:,:,kdir),b_k_t(:,:,kdir),b_k(:,:,jdir),b_k_t(:,:,jdir))

                b_k( :,:,kdir) = b_k( :,:,kdir) - mu * b_k( :,:,jdir)
                dQ_k(:,:,kdir) = dQ_k(:,:,kdir) - mu * dQ_k(:,:,jdir)
                if (iflow_type > FLOW_LAMINAR) then
                    b_k_t( :,:,kdir) = b_k_t( :,:,kdir) - mu * b_k_t( :,:,jdir)
                    dQ_k_t(:,:,kdir) = dQ_k_t(:,:,kdir) - mu * dQ_k_t(:,:,jdir)
                endif
                norm_b_k_inv = one / l2norm(nq,nturb,ncells,b_k(:,:,kdir),b_k_t(:,:,kdir))
                
                b_k( :,:,kdir) = b_k( :,:,kdir) * norm_b_k_inv
                dQ_k(:,:,kdir) = dQ_k(:,:,kdir) * norm_b_k_inv
                if (iflow_type > FLOW_LAMINAR) then
                    dQ_k(  :,:,kdir) = dQ_k(  :,:,kdir) * norm_b_k_inv
                    dQ_k_t(:,:,kdir) = dQ_k_t(:,:,kdir) * norm_b_k_inv
                endif
            enddo

            ! Update correction and residual
            gamma_k = inner_product(nq,nturb,ncells,b_k(:,:,kdir),b_k_t(:,:,kdir),r_k,r_k_t) ! r_k is still r_(k-1) at this point
            
            gcr_final_update_f = gcr_final_update_f + gamma_k * dQ_k(:,:,kdir)
            
            r_k = r_k - gamma_k * b_k(:,:,kdir) ! r_k is now up to date
            if (iflow_type > FLOW_LAMINAR) then
                gcr_final_update_T = gcr_final_update_T + gamma_k * dQ_k_t(:,:,kdir)
                r_k_t = r_k_t - gamma_k * b_k_t(:,:,kdir) ! r_k_t is now up to date
            endif
            ! Check for convergence
            rms_r_k = rms(nq,nturb,ncells,r_k,r_k_t,inv_ncells)
            if ( ( rms_r_k / rms_r_0 ) < gcr_reduction_target ) then
                iostat = GCR_SUCCESS
                n_projections = kdir
                nl_reduction  = rms_r_k / rms_r_0
                gcr_res_rms   = rms_r_k
                ! Clear Jacobian arrays
                call clear_jacobian_arrays(V,C,R,Dinv)
                call clear_jacobian_arrays(Vt,Ct,Rt,Dinvt)
                return
            end if

            ! Check for stall
            ! Original method
            norm_r_k = l2norm(nq,nturb,ncells,r_k,r_k_t)
            if (gamma_k < norm_r_k * 0.001_p2) then
                iostat = GCR_STALL
                n_projections = jdir
                ! Clear Jacobian arrays
                call clear_jacobian_arrays(V,C,R,Dinv)
                call clear_jacobian_arrays(Vt,Ct,Rt,Dinvt)
                return
            endif
        end do proj_loop

        ! if we make it this far we've stalled
        iostat = GCR_STALL
        n_projections = jdir
        ! Clear Jacobian arrays
        call clear_jacobian_arrays(V,C,R,Dinv)
        call clear_jacobian_arrays(Vt,Ct,Rt,Dinvt)
        return
    end subroutine gcr_solve_scratch

    subroutine compute_frechet(q,res,qT,resT,dQ_f,dQ_t,mag_dQ,sol_rms,frechet_deriv_f,frechet_deriv_t,os)

        use common , only : p2, one, half

        use utils , only : iflow_type, FLOW_LAMINAR

        use solution_vars , only : nq, dtau

        use turb , only : nturb
        
        use solution , only : compute_primative_jacobian

        use config , only : gcr_verbosity, CFL_turb

        use grid , only : ncells, cell

        use residual , only : compute_residual

        use turb , only : twsn

        implicit none

        real(p2), dimension(:,:), pointer, intent(inout) :: q, res ! Flow values
        real(p2), dimension(:,:), pointer, intent(inout) :: qT, resT ! Turbulent Flow values
        real(p2), dimension(:,:),          intent(in)    :: dQ_f
        real(p2), dimension(:,:),          intent(in)    :: dQ_t
        real(p2),                          intent(in)    :: mag_dQ
        real(p2),                          intent(in)    :: sol_rms
        real(p2), dimension(:,:),          intent(out)   :: frechet_deriv_f
        real(p2), dimension(:,:),          intent(out)   :: frechet_deriv_t
        integer,                           intent(out)   :: os

        real(p2), dimension(:,:), pointer :: q_n, r_0
        real(p2), dimension(:,:), pointer :: t_n, tr_0
        real(p2), dimension(5,5)          :: prim_jac
        real(p2)                          :: eps_frechet
        real(p2)                          :: frech_min_bound = 1.0e-07_p2
        real(p2), dimension(2)            :: dtaui

        integer :: icell, it

        os = 0

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

        do icell = 1,ncells
            q(:,icell) = q_n(:,icell) + eps_frechet * dQ_f(:,icell) / mag_dQ
            if ( (q(1,icell) <= 0) .or. ( q(5,icell) <= 0 ) ) then
                ! solution does not pass realizeability check
                deallocate(  q)
                deallocate(res)

                q   => q_n
                res => r_0
                
                os = GCR_PREAL_FAIL
                return
            endif
        end do
        if (iflow_type > FLOW_LAMINAR) then
            ! move the lates solution vector to the temp vector
            t_n  => qT
            tr_0 => resT
            nullify(qT, resT)
            allocate(qT(ncells,nturb))
            allocate(resT(ncells,nturb))
            ltrb: do it = 1,nturb 
                do icell = 1,ncells
                    qT(icell,it) = t_n(icell,it) + eps_frechet * dQ_t(icell,it) / mag_dQ
                end do 
            end do ltrb
            ! Currently there are no realizeability checks for the turbulent variables
        end if

        call compute_residual

        frechet_deriv_f = mag_dQ * ( res - r_0 ) / eps_frechet
        
        ! write(*,*) "  Pre psuedo time:", l2norm(nq,ncells,frechet_deriv(:,:))
        ! Nevermind its equation 14
        ! it including it seems incorrect
        do icell = 1,ncells
            prim_jac = compute_primative_jacobian(q(:,icell))
            frechet_deriv_f(:,icell) = frechet_deriv_f(:,icell) + cell(icell)%vol/dtau(icell) * matmul(prim_jac,dQ_f(:,icell))
            ! frechet_deriv(:,icell) = frechet_deriv(:,icell) + cell(icell)%vol/dtau(icell) * sol_update(:,icell)
        end do
        ! write(*,*) " Post psuedo time:", l2norm(nq,ncells,frechet_deriv(:,:))

        deallocate(  q)
        deallocate(res)

        q   => q_n
        res => r_0

        ! nullify(q_n,r_0) ! These are nullified when they go out of scope

        if (iflow_type > FLOW_LAMINAR) then 
            frechet_deriv_t = mag_dQ * ( resT - tr_0 ) / eps_frechet
            do it = 1,nturb
                do icell = 1,ncells
                    dtaui(1) = CFL_turb * cell(icell)%vol/( half * twsn(1,icell) )
                    dtaui(2) = CFL_turb * (cell(icell)%vol)**2 / (twsn(2,icell))
                    qT(icell,it) = qT(icell,it) + dQ_t(icell,it) * cell(icell)%vol / minval(dtaui)
                end do
            end do
            deallocate(resT,qT)
            qT => t_n
            resT => tr_0
            ! nullify(t_n, tr_0)
        end if

    end subroutine compute_frechet

    subroutine gcr_real_check(sol_current,sol_update, iostat)

        ! Realizability chekc for the proposed nonlinear update from gcr_run

        use common  , only : p2

        use solution_vars, only : nq

        use grid    , only : ncells

        implicit none 

        real(p2), dimension(nq,ncells), intent(in) :: sol_current
        real(p2), dimension(nq,ncells), intent(in) :: sol_update
        integer                       , intent(out):: iostat

        integer :: icell
        
        do icell = 1,ncells
            if ( .not.check_non_real_update(nq, sol_current(:,icell), sol_update(:,icell)) ) then
                iostat = GCR_REAL_FAIL
                return
            endif
        end do  

        ! If we made it this far we know the update was valid
        iostat = GCR_SUCCESS

    end subroutine gcr_real_check

    pure function check_non_real_update(nq,sol_current,sol_update) result(is_real)
        !
        ! Check that the proposed solution does not result in a non-real solution (such as negative pressure or temperature) 
        ! Output:
        ! .TRUE. = The proposed solution real (pass)
        ! .FALSE.= The proposed solution is not real (fail)
        use common , only : p2, zero

        implicit none

        integer ,                intent(in) :: nq
        real(p2), dimension(nq), intent(in) :: sol_current
        real(p2), dimension(nq), intent(in) :: sol_update
        logical                             :: is_real

        ! if ( (sol_current(1) + sol_update(1) <= zero ) .OR. (sol_current(5) <= sol_update(5)) ) then
        !     is_real = .false.
        ! else
        !     is_real = .true.
        ! end if
        
        is_real = ( (sol_current(1) + sol_update(1) > zero ) .AND. (sol_current(5) + sol_update(5) > zero) )
    end function check_non_real_update

    subroutine gcr_nl_control(sol_update_f,sol_update_t, gcr_res_rms, iostat)

        ! Determine the optimum underrelaxation factor for the nonlinear solution update.

        use common , only : p2, half, one, two

        use utils , only : iflow_type, FLOW_LAMINAR

        use config , only : gcr_reduction_target, CFL_turb

        use solution_vars , only : nq, q, res, dtau, inv_ncells

        use solution , only : compute_primative_jacobian

        use grid , only : ncells, cell

        use residual , only : compute_residual

        use turb , only : nturb, turb_var, turb_res, twsn

        implicit none

        real(p2), dimension(nq,ncells),    intent(in) :: sol_update_f
        real(p2), dimension(ncells,nturb), intent(in) :: sol_update_t
        real(p2),                          intent(in) :: gcr_res_rms
        integer,                           intent(out):: iostat

        real(p2), dimension(:,:), pointer   :: q_n, r_0
        real(p2), dimension(:,:), pointer   :: t_n, tr_0
        real(p2)                            :: residual_reduct_target
        real(p2)                            :: Rtau_rms, R0_rms
        real(p2)                            :: delQ_rms, Qn_rms
        real(p2)                            :: f_0, f_1, g_1     ! terms in the optimization equation 22
        real(p2), dimension(nq,ncells)      :: frechet_deriv_f
        real(p2), dimension(ncells,nturb)   :: frechet_deriv_t
        real(p2)                            :: delQ_norm
        real(p2)                            :: ur_opt, ur_min
        real(p2), dimension(2)              :: dtaui

        integer :: icell, it

        q_n => q
        r_0 => res

        nullify(q,res)

        allocate(  q(nq,ncells))
        allocate(res(nq,ncells))

        ! We don't need to perform a realizability check here because we did it in gcr_real_check
        q = q_n + sol_update_f

        if (iflow_type > FLOW_LAMINAR) then
            t_n  => turb_var
            tr_0 => turb_res
            nullify(turb_var, turb_res)
            allocate(turb_var(ncells,nturb))
            allocate(turb_res(ncells,nturb))
            turb_var = t_n + sol_update_t
        endif
        

        call compute_residual

        do icell = 1,ncells
            res(:,icell) = res(:,icell) + cell(icell)%vol/dtau(icell) *  &
                            matmul( compute_primative_jacobian(q(:,icell)) , sol_update_f(:,icell) )
        end do  

        do it = 1,nturb
            do icell = 1,ncells

                ! TODO add pseudo-transient term to this
                dtaui(1) = CFL_turb * cell(icell)%vol/( half * twsn(1,icell) )
                dtaui(2) = CFL_turb * (cell(icell)%vol)**2 / (twsn(2,icell))
                turb_res(icell,it) = turb_res(icell,it) + sol_update_t(icell,it) * cell(icell)%vol / minval(dtaui)
            end do
        end do


        residual_reduct_target = half * (one + gcr_reduction_target)

        Rtau_rms = rms(nq,nturb,ncells,res,turb_res,inv_ncells)
        R0_rms   = rms(nq,nturb,ncells,r_0,tr_0,inv_ncells)

        if ( Rtau_rms / R0_rms < residual_reduct_target ) then
            ! The reduction is acceptable without underrelaxation
            ! Additionally q and res have already been updated.
            deallocate(q_n)
            deallocate(r_0)
            if (iflow_type > FLOW_LAMINAR) then
                deallocate(t_n)
                deallocate(tr_0)
            endif
            iostat = GCR_SUCCESS
            return
        endif

        ! If the change isn't succesful first we will check if the change is comperable to the computer percision
        delQ_rms = rms(nq,nturb,ncells,sol_update_f,sol_update_t,inv_ncells)
        Qn_rms   = rms(nq,nturb,ncells,q_n         ,t_n         ,inv_ncells)

        ! write(*,"(a,es12.6)") "delQ_rms/Qn_rms = ", delQ_rms / Qn_rms 
        if ( delQ_rms / Qn_rms < 1.0e-12_p2) then
            if ( Rtau_rms / R0_rms < one) then
                ! Any reduction will be considered a success at this point, but we don't want to make the CFL any bigger
                deallocate(q_n)
                deallocate(r_0)
                if (iflow_type > FLOW_LAMINAR) then
                    deallocate(t_n)
                    deallocate(tr_0)
                endif
                iostat = GCR_CFL_FREEZE
                return
            endif
        endif

        ! The residual did not reduce so now we will apply an underrelaxation factor to minimize the residual
        f_0 = R0_rms
        f_1 = Rtau_rms

        ! At this point:
        ! q_n, r_0 => Previous (outer) iteration value
        ! q, res   => Proposed updated solution (w/o underrelaxation)
        ! ditto for the turb variables

        ! We want the frechet derivitive evaluated at the previous solution
        delQ_norm = l2norm(nq,nturb,ncells,sol_update_f, sol_update_t)
        call compute_frechet(q_n,r_0,t_n,tr_0,sol_update_f,sol_update_t,delQ_norm,Qn_rms,frechet_deriv_f,frechet_deriv_t,iostat)

        ! TODO: add error handling for nonzero frechet iostat
        
        ! We will temporarily reuse the res vector to save memory space
        do icell = 1,ncells
            ! EQ 21 from FUN 3D paper where omega = 1
            res(:,icell) = r_0(:,icell) + matmul( compute_primative_jacobian(q(:,icell)) , sol_update_f(:,icell) ) *&
                            cell(icell)%vol/dtau(icell) + frechet_deriv_f(:,icell)
        end do
        do it = 1,nturb
            do icell = 1,ncells
                dtaui(1) = CFL_turb * cell(icell)%vol/( half * twsn(1,icell) )
                dtaui(2) = CFL_turb * (cell(icell)%vol)**2 / (twsn(2,icell))
                turb_res(icell,it) = tr_0(icell,it) + cell(icell)%vol / minval(dtaui) * sol_update_t(icell,it) &
                                     + frechet_deriv_t(icell,it)
            end do
        end do
        g_1 = rms(nq,nturb,ncells,res,turb_res,inv_ncells)
        
        ! We already have the g_1 term from our last projection
        ! g_1 = gcr_res_rms ! No we don't...

        ! Minimize the quadratic function a*ur**2 + b*ur + c
        ! where c = f_0, b = g_1-c, and a = f_1-b-c.
        ur_opt = ( g_1 - f_0 ) / (two * (f_1 - g_1) )
        ! Adding a bound to the ur factor 0 <= ur <= 1
        ur_min = ( one-residual_reduct_target ) / ( one - (gcr_res_rms/R0_rms) )
        ur_opt = min(max(ur_opt,ur_min) , one)

        ! Because this passed the realizability check with ur = 1, we know -Q(j,i) < sol_update for j = 1,5 and any i.
        ! Therefore if abs(ur_opt) < 1, we now the updated solution w/ under-relaxation will also be realizable.
        q = q_n + ur_opt * sol_update_f
        if (iflow_type > FLOW_LAMINAR) turb_var = t_n + ur_opt * sol_update_t

        ! Check convergence of the updated solution
        call compute_residual

        ! Compute R_tau
        do icell = 1,ncells
            res(:,icell) = res(:,icell) + cell(icell)%vol/dtau(icell) * &
                        matmul( compute_primative_jacobian(q(:,icell)) , sol_update_f(:,icell) )
        end do  
        do it = 1,nturb
            do icell = 1,ncells
                dtaui(1) = CFL_turb * cell(icell)%vol/( half * twsn(1,icell) )
                dtaui(2) = CFL_turb * (cell(icell)%vol)**2 / (twsn(2,icell))
                turb_res(icell,it) = turb_res(icell,it) + cell(icell)%vol / minval(dtaui) * sol_update_t(icell,it)
            end do
        end do
        Rtau_rms = rms(nq,nturb,ncells,res,turb_res,inv_ncells)

        if ( Rtau_rms / R0_rms < residual_reduct_target .OR. delQ_rms / Qn_rms < 1.0e-12_p2) then
            ! q and res have already been updated and the residual has reduced.
            iostat = GCR_CFL_FREEZE
            ! along some of the paths the qcombs get pointed back to their counter parts.  In this case we need to just nullify
            deallocate(q_n)
            deallocate(r_0)
            if (iflow_type > FLOW_LAMINAR) then
                deallocate(t_n)
                deallocate(tr_0)
            endif
            return
        endif

        ! If we've made it this far we failed :(
        ! Undo the solution update and report failure
        ! Revert the solution and residual to retry
        deallocate(q)
        deallocate(res)
        q   => q_n ! q_n will be nullified when it goes out of scope
        res => r_0 ! r_0 will be nullified when it goes out of scope
            if (iflow_type > FLOW_LAMINAR) then
                deallocate(turb_res)
                deallocate(turb_var)
                turb_var   => t_n
                turb_res => tr_0
            endif
        iostat = GCR_STALL

    end subroutine gcr_nl_control

    subroutine gcr_CFL_control(gcr_status)

        use common          , only : p2, two, half
        
        use config          , only : CFL, CFL_max, CFL_min, CFL_turb

        use grid            , only : ncells, cell

        use solution_vars   , only : jac, q, nq, dtau, CFL_used

        use solution        , only : compute_primative_jacobian, compute_local_time_step_dtau

        use direct_solve    , only : gewp_solve, safe_invert_scalar

        use turb            , only : nturb, turb_jac, twsn, turb_var

        implicit none

        integer, intent(inout) :: gcr_status ! not actually used for now...

        real(p2), dimension(nq,nq) :: preconditioner

        integer :: icell, k, j, it
        integer :: idestat

        real(p2), dimension(2) :: dtaui


        if (gcr_status == GCR_SUCCESS) then
            if (gcr_verbosity >= 3) then
                write(*,*) "Successful iteration:"
                write(*,*) "CFL used:   ", CFL
            endif
            CFL_used = CFL
            CFL = min(CFL * two, CFL_max)
            CFL_turb = CFL
            if (gcr_verbosity >= 3) then
                write(*,*) "CFL update: ", CFL
            endif
        elseif (gcr_status == GCR_CFL_FREEZE) then
            ! Nothing to actually do here other than report a successful update
            CFL_used = CFL
            gcr_status = GCR_SUCCESS
        else ! fail
            do icell = 1,ncells
                preconditioner = compute_primative_jacobian(q(:,icell))

                ! we want to remove the pseudo transient term so that we can add it with a different CFL
                jac(icell)%diag = jac(icell)%diag - (cell(icell)%vol/dtau(icell))*preconditioner
                
            end do

            do it = 1,nturb
                do icell = 1,ncells
                    dtaui(1) = CFL_turb * cell(icell)%vol/( half * twsn(1,icell) )
                    dtaui(2) = CFL_turb * (cell(icell)%vol)**2 / (twsn(2,icell))
                    turb_jac(icell,it)%diag = turb_jac(icell,it)%diag - cell(icell)%vol / minval(dtaui) * turb_var(icell,it)
                end do
            end do

            if (gcr_verbosity >= 3) then
                write(*,*) "Stall detected:"
                write(*,*) "CFL old: ", CFL
            endif
            CFL = max(CFL / 10.0_p2, CFL_min)
            CFL_turb = CFL
            if (gcr_verbosity >= 3) then
                write(*,*) "CFL new: ", CFL
            endif
            ! Update the time step with the new CFL
            call compute_local_time_step_dtau

            ! Update the jac_diag and jac_diag_inv
            do icell = 1,ncells
                preconditioner = compute_primative_jacobian(q(:,icell))

                ! we want to remove the pseudo transient term so that we can add it with a different CFL
                jac(icell)%diag = jac(icell)%diag + (cell(icell)%vol/dtau(icell))*preconditioner
                
                ! Invert the diagonal
                idestat = 0
                !                A                     dim  A^{-1}               error check
                call gewp_solve( jac(icell)%diag(:,:), 5  , jac(icell)%diag_inv, idestat    )
                !  Report errors
                if (idestat/=0) then
                    write(*,*) " Error in inverting the diagonal block... Stop"
                    write(*,*) "  Cell number = ", icell
                    do k = 1, 5
                        write(*,'(12(es8.1))') ( jac(icell)%diag(k,j), j=1,5 )
                    end do
                    stop
                endif

            end do

            do it = 1,nturb
                do icell = 1,ncells
                    dtaui(1) = CFL_turb * cell(icell)%vol/( half * twsn(1,icell) )
                    dtaui(2) = CFL_turb * (cell(icell)%vol)**2 / (twsn(2,icell))
                    turb_jac(icell,it)%diag = turb_jac(icell,it)%diag + cell(icell)%vol / minval(dtaui) * turb_var(icell,it)
                    turb_jac(icell,1)%diag_inv = safe_invert_scalar(turb_jac(icell,1)%diag)
                end do
            end do
        endif

    end subroutine gcr_CFL_control

    pure function rms(nq,nturb,ncells,Vflow, Vturb,div)

        ! Function for computing the L2 norm

        use common , only : p2,zero 

        implicit none

        integer, intent(in)                 :: nq, nturb, ncells
        real(p2),dimension(:,:), intent(in) :: Vflow
        real(p2),dimension(:,:), intent(in) :: Vturb
        real(p2),                intent(in) :: div     ! dividend
        real(p2)                            :: rms

        integer :: i,j

        rms = zero

        ! Compute flow variables
        do i = 1,ncells
            rms = rms + dot_product(Vflow(:,i),Vflow(:,i))
        end do

        do i = 1,nturb
            rms = rms + dot_product(Vturb(:,i),Vturb(:,i)) ! this works lmao
        end do

        rms = sqrt(rms*div)
    end function

    pure function l2norm(nq,nturb,ncells,Vflow, Vturb)

        ! Function for computing the L2 norm of block vectors

        use common , only : p2,zero 

        implicit none

        integer, intent(in)                 :: nq, nturb, ncells
        real(p2),dimension(:,:), intent(in) :: Vflow
        real(p2),dimension(:,:), intent(in) :: Vturb
        real(p2)                            :: l2norm

        integer :: i,j

        l2norm = zero

        do i = 1,ncells
            l2norm = l2norm + dot_product(Vflow(:,i),Vflow(:,i))
        end do

        do i = 1,nturb
            l2norm = l2norm + dot_product(Vturb(:,i),Vturb(:,i)) ! this works lmao
        end do

        l2norm = sqrt(l2norm)
    end function

    pure function inner_product(nq,nturb,ncells,vector1F,vector1T,vector2F,vector2T)

        ! Function for computing the inner_product of block vectors

        use common , only : p2,zero 

        implicit none

        integer, intent(in)                 :: nq, nturb, ncells
        real(p2),dimension(:,:), intent(in) :: vector1F, vector2F
        real(p2),dimension(:,:), intent(in) :: vector1T, vector2T
        real(p2)                            :: inner_product

        integer :: i,j

        inner_product = zero


        do i = 1,ncells
            inner_product = inner_product + dot_product( vector1F(:,i) , vector2F(:,i) )
        end do
        do i = 1,nturb
            inner_product = inner_product + dot_product( vector1T(:,i) , vector2T(:,i) )
        end do

    end function

    subroutine clear_jacobian_arrays_B(V,C,R,Dinv)

        use common , only : p2

        implicit none

        real(p2), dimension(:,:,:), pointer, intent(inout) :: V   ! Values (5x5 block matrix) plus corresponding index
        integer , dimension(:),     pointer, intent(inout) :: C   ! Column index of each value
        integer , dimension(:),     pointer, intent(inout) :: R   ! Start index of each new row
        real(p2), dimension(:,:,:), pointer, intent(inout) :: Dinv

        if (associated(V)) deallocate(V)
        if (associated(C)) deallocate(C)
        if (associated(R)) deallocate(R)
        if (associated(Dinv)) deallocate(Dinv)

    end subroutine clear_jacobian_arrays_B

    subroutine clear_jacobian_arrays_S(V,C,R,Dinv)

        use common , only : p2

        implicit none

        real(p2), dimension(:,:), pointer, intent(inout) :: V   ! Values (5x5 block matrix) plus corresponding index
        integer , dimension(:),     pointer, intent(inout) :: C   ! Column index of each value
        integer , dimension(:),     pointer, intent(inout) :: R   ! Start index of each new row
        real(p2), dimension(:,:), pointer, intent(inout) :: Dinv

        if (associated(V)) deallocate(V)
        if (associated(C)) deallocate(C)
        if (associated(R)) deallocate(R)
        if (associated(Dinv)) deallocate(Dinv)

    end subroutine clear_jacobian_arrays_S
    
end module gcr