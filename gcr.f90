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



    contains

    subroutine gcr_run(iostat)

        use common      , only : p2

        use grid        , only : ncells

        use solution    , only : nq, jacobian_type, q

        implicit none
        
        integer ,          intent(out) :: iostat

        real(p2), dimension(nq,ncells) :: sol_update
        real(p2)                       :: gcr_res_rms

        call gcr_solve(sol_update,gcr_res_rms,iostat)

        if (iostat /= GCR_SUCCESS) return

        call gcr_real_check(q,sol_update,iostat)

        if (iostat /= GCR_SUCCESS) return

        call gcr_nl_control(sol_update,gcr_res_rms,iostat)

    end subroutine gcr_run

    subroutine gcr_solve(gcr_final_update, gcr_res_rms, iostat)

        ! Right preconditioned version of Algorithm 6.21 from Y. Saad. Iterative Methods for Sparse Linear Systems, 2nd Edition
        ! Right preconditioning is done using the same linear solver and the approximate (1st order jacobian)

        use common      , only : p2, zero, my_eps, one

        use config      , only : gcr_max_projections, gcr_reduction_target, amg_cycle

        use grid        , only : ncells, cell

        use solution    , only : nq, inv_ncells, compute_primative_jacobian, jacobian_type, res, jac, &
                                 nl_reduction, n_projections, q, roc

        use residual    , only : compute_residual

        use linear_solver, only: build_A_BCSM, build_Dinv_array, linear_sweeps, RELAX_FAIL_STALL, RELAX_FAIL_DIVERGE, & 
                                 multilevel_cycle

        use algebraic_multigird, only : convert_amg_c_to_i

        implicit none

        real(p2),            dimension(nq,ncells), intent(out) :: gcr_final_update   ! delta_Q_n+1
        real(p2),                                  intent(out) :: gcr_res_rms        ! Needed for the NL update subroutine
        integer,                                   intent(out) :: iostat                 ! status of gcr solve

        !+++++++++++++++++
        real(p2), dimension(nq,ncells)                     :: gcr_residual 
        real(p2)                                           :: RMS_R0, RMS_Qn ! RMS of outer residual and solution vectors
        real(p2)                                           :: rms_resj, rms_resj_m1
        real(p2), dimension(nq,ncells,gcr_max_projections) :: p ! Preconditioned vector p_i = M^{-1}*r_i
        real(p2), dimension(nq,ncells,gcr_max_projections) :: Ap ! Exact jacobian (A computed w Frechet Deriv) times p 
        real(p2)                                           :: length_p
        real(p2), dimension(gcr_max_projections)           :: length_Ap2
        real(p2)                                           :: alpha, beta
        logical                                            :: stall_cond
        real(p2)                                           :: gammak_maxApk
        real(p2)                                           :: single_reduct

        ! Variables for preconditioning matrix M
        real(p2), dimension(:,:,:), pointer     :: V   ! Values (5x5 block matrix) plus corresponding index
        integer , dimension(:),     pointer     :: C   ! Column index of each value
        integer , dimension(:),     pointer     :: R   ! Start index of each new row
        integer                                 :: nnz
        real(p2), dimension(:,:,:), pointer     :: Dinv
        integer                                 :: direction, level

        integer :: idir, jdir
        integer :: cycle_type
        integer :: os

        gcr_final_update = zero
        gcr_residual = -res
        gcr_res_rms = zero
        RMS_R0 = rms(nq,ncells,res,inv_ncells)
        RMS_Qn = rms(nq,ncells,q  ,inv_ncells)
        stall_cond = .false.
        rms_resj = zero
        single_reduct = gcr_reduction_target ** (one/real(gcr_max_projections,p2))

        jdir = 1

        if (gcr_verbosity >= 3) write(*,"(a,es12.6)") "RMS_R0 = ", RMS_R0 

        ! Build M (A Approx) for precondition solve
        allocate(R(ncells+1))
        allocate(Dinv(5,5,ncells))
        call build_A_BCSM(ncells,cell,jac,V,C,R,nnz=nnz)

        call build_Dinv_array(ncells,jac,Dinv)

        level = 1

        ! Compute the initial search direction
        cycle_type = convert_amg_c_to_i(amg_cycle)
        ! keep_A = .true. so that V,C,R, and Dinv do not have to be rebuilt
        call multilevel_cycle(ncells,nq, V, C, R, -gcr_residual, Dinv,cycle_type,.true.,p(:,:,jdir),os)
        ! call linear_sweeps(ncells,nq,nnz,V,C,R,-gcr_residual,Dinv,level,direction,p(:,:,jdir),os)
        if (os == RELAX_FAIL_DIVERGE) then
            iostat = GCR_PRECOND_DIVERGE
            return
        elseif (os == RELAX_FAIL_STALL) then
            iostat = GCR_PRECOND_STALL
            return
        endif
        length_p = l2norm(nq,ncells,p(:,:,jdir))
        
        call compute_frechet(p(:,:,jdir),length_p,RMS_Qn,Ap(:,:,jdir))


        proj_loop : do

            length_Ap2(jdir) = max(l2norm(nq,ncells,Ap(:,:,jdir))**2,my_eps)

            alpha = inner_product(nq,ncells,gcr_residual(:,:),Ap(:,:,jdir)) / (length_Ap2(jdir))

            gcr_final_update = gcr_final_update + alpha * p(:,:,jdir)

            gcr_residual     = gcr_residual     - alpha * Ap(:,:,jdir)

            ! Check for sufficient rms reduction
            rms_resj_m1 = rms_resj
            rms_resj = rms(nq,ncells,gcr_residual,inv_ncells)
            if (rms_resj / RMS_R0 < gcr_reduction_target) then
                iostat = GCR_SUCCESS
                n_projections = jdir
                nl_reduction  = rms_resj / RMS_R0
                gcr_res_rms   = rms_resj
                return
            ! Check for max projections
            elseif (jdir >= gcr_max_projections) then
                iostat = GCR_STALL
                n_projections = jdir
                return
            endif
            ! Check for stall
            ! https://doi.org/10.2514/6.2019-2333 Eq (7)
            ! gammak = (Ap^{k} , r^{k-1}) = alpha * | Ap^{k} |**2
            gammak_maxApk = maxval(Ap(:,:,jdir))
            gammak_maxApk = gammak_maxApk * alpha * length_Ap2(jdir)
            if (jdir > 1) then
                if (gammak_maxApk < rms_resj) then
                    ! Potential stall
                    if (rms_resj/rms_resj_m1 > single_reduct) then
                        iostat = GCR_STALL
                        n_projections = jdir
                        return
                    endif
                endif
            endif
            ! ! Check for stall
            ! if (jdir > 1) then
            !     if (gcr_verbosity >= 3) then
            !         write(*,*) "res * Ap_jdir", inner_product(nq,ncells,gcr_residual(:,:),p(:,:,idir))
            !     endif
            !     stall_cond = alpha * maxval(Ap(:,:,jdir)) > rms_resj
            ! endif
            ! if (stall_cond) then
            !     iostat = GCR_STALL
            !     n_projections = jdir
            !     return
            ! endif


            ! Generate new search direction
            call multilevel_cycle(ncells,nq, V, C, R, -gcr_residual, Dinv,cycle_type,.true.,p(:,:,jdir+1),os)
            ! call linear_sweeps(ncells,nq,nnz,V,C,R,-gcr_residual,Dinv,level,direction,p(:,:,jdir+1),os)
            if (os == RELAX_FAIL_DIVERGE) then
                iostat = GCR_PRECOND_DIVERGE
                return
            elseif (os == RELAX_FAIL_STALL) then
                iostat = GCR_PRECOND_STALL
                return
            endif
            
            call compute_frechet(p(:,:,jdir+1),length_p,RMS_Qn,Ap(:,:,jdir+1))

            ! Modified Gram-Schmidt Orthogonalization
            do idir = 1,jdir
                beta = - inner_product(nq,ncells,Ap(:,:,jdir+1),Ap(:,:,idir)) / length_Ap2(idir)
                
                p(:,:,jdir+1) = p(:,:,jdir+1) + beta * p(:,:,idir)
                Ap(:,:,jdir+1)=Ap(:,:,jdir+1) + beta *Ap(:,:,idir)

                ! Sanity check orthogonality
                if (gcr_verbosity >= 3) then
                    ! Note: the P vectors are A^{T}A-orthogonal. meaning Ap_i * Ap_j = 0.  However, the vectors P_i * P_j need not be 
                    ! orthogonal themselves.  In fact, they won't be unless A^{T}A = c*I (I think, I haven't done an exhaustive 
                    ! proof)
                    if (idir == 1) write(*,"(a,i2,a,es12.6)") " Res(", jdir,") = ", rms_resj
                    write(*,*) "j=",jdir+1,"i=",idir,"inner_product (p) = ", inner_product(nq,ncells,p(:,:,jdir + 1),p(:,:,idir))
                    write(*,*) "j=",jdir+1,"i=",idir,"inner_product (Ap) = ", inner_product(nq,ncells,Ap(:,:,jdir + 1),Ap(:,:,idir))
                endif
                continue
            end do

            jdir = jdir + 1
                
        end do proj_loop

    end subroutine gcr_solve

    subroutine gcr_solve_scratch(gcr_final_update, iostat)

        use common      , only : p2, zero

        use config      , only : gcr_max_projections, gcr_reduction_target

        use grid        , only : ncells, cell

        use solution    , only : nq, inv_ncells, dtau, compute_primative_jacobian, jacobian_type, res, jac, &
                                 nl_reduction, n_projections, q

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
        real(p2)                                           :: gamma_k ! gamma_k from EQ. 14
        real(p2), dimension(nq,ncells)                     :: gcr_residual
        real(p2)                                           :: gcr_res_rms
        real(p2)                                           :: max_bk



        integer :: kdir, jdir ! projection direction indices
        integer :: icell
        integer :: os


        ! Debug:
        if (gcr_verbosity >= 4) write(*,*) "Q(:,5): ", q(:,5)
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
                ! I originally included the primitive jacobian (dW/dQ) in the psuedo-time term but I don't think I'm supposed to 
                ! include it for two reasons:
                ! 1) it doesn't converge if I include it
                ! 2) delta_Q is already the update to Q as opposed to an update to W, so we shouldn't be converting the update

                b_k(:,icell,kdir) = b_k(:,icell,kdir) + delta_Q_k(:,icell,kdir) * cell(icell)%vol/dtau(icell)
            end do

            b_k_mag = l2norm(nq,ncells,b_k(:,:,kdir))

            ! Normalize the correction and search direction
            delta_Q_k(:,:,kdir) = delta_Q_k(:,:,kdir) / b_k_mag

            b_k(:,:,kdir)       = b_k(:,:,kdir)       / b_k_mag

            ! Normalize direction k to previous search direactions
            jdir = 1

            do while (jdir < kdir)
                ! Pretty sure this is just a modified Gram-Schmidt... suppose I could check...
                gcr_inner_prod = inner_product(nq,ncells, b_k(:,:,kdir), b_k(:,:,jdir)) ! mu

                b_k(:,:,kdir)       = b_k(:,:,kdir)       - gcr_inner_prod * b_k(:,:,jdir)
                delta_Q_k(:,:,kdir) = delta_Q_k(:,:,kdir) - gcr_inner_prod * delta_Q_k(:,:,jdir)

                ! Renormalize
                b_k_mag = l2norm(nq,ncells,b_k(:,:,kdir)) 

                delta_Q_k(:,:,kdir) = delta_Q_k(:,:,kdir) / b_k_mag
                b_k(:,:,kdir)       = b_k(:,:,kdir)       / b_k_mag

                jdir = jdir + 1
            end do

            ! Compute projection of the current search direction and the previous residual
            gamma_k = inner_product(nq,ncells, b_k(:,:,kdir), gcr_residual(:,:))

            ! if they are (nearly) orthogonal the solution has stalled and we should quit
            gcr_res_rms = rms(nq,ncells,gcr_residual,inv_ncells)
            max_bk = maxval(b_k(:,:,kdir))
            ! if (gcr_verbosity >= 3 ) write(*,*) "max_bk: ", max_bk
            max_bk = gamma_k * max_bk
            ! if (gcr_verbosity >= 3 ) then
            !     write(*,*) "gamma_k * max_bk: ", max_bk
            !     write(*,*) "gcr_res_rms: ", gcr_res_rms
            ! endif
            
            ! if ( max_bk < gcr_res_rms ) exit project_loop

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
        if (gcr_verbosity >= 3) write(*,*) "STALL: Reduction: ", nl_reduction

    end subroutine gcr_solve_scratch

    subroutine compute_frechet(sol_update,update_length,sol_rms,frechet_deriv)

        use common , only : p2, one

        use solution , only : nq, q, res, dtau, compute_primative_jacobian

        use grid , only : ncells, cell

        use residual , only : compute_residual

        implicit none

        real(p2), dimension(nq,ncells), intent(in) :: sol_update
        real(p2),                       intent(in) :: update_length
        real(p2),                       intent(in) :: sol_rms
        real(p2), dimension(nq,ncells), intent(out):: frechet_deriv

        real(p2), dimension(:,:), pointer :: q_n
        real(p2), dimension(:,:), pointer :: r_0
        real(p2), dimension(5,5)          :: prim_jac
        real(p2)                          :: eps_frechet
        real(p2)                          :: frech_min_bound = 1.0e-07_p2

        integer :: icell

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

        do icell = 1,ncells
            prim_jac = compute_primative_jacobian(q_n(:,icell))
            frechet_deriv(:,icell) = frechet_deriv(:,icell) + cell(icell)%vol/dtau(icell) * matmul(prim_jac,sol_update(:,icell))
            ! frechet_deriv(:,icell) = frechet_deriv(:,icell) + cell(icell)%vol/dtau(icell) * sol_update(:,icell)
        end do

        deallocate(  q)
        deallocate(res)

        q   => q_n
        res => r_0

        nullify(q_n,r_0)

    end subroutine compute_frechet

    subroutine gcr_real_check(sol_current,sol_update, iostat)

        ! Realizability chekc for the proposed nonlinear update from gcr_run

        use common  , only : p2

        use solution, only : nq

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

    subroutine gcr_nl_control(sol_update, gcr_res_rms, iostat)

        ! Determine the optimum underrelaxation factor for the nonlinear solution update.

        use common , only : p2, half, one, two, zero

        use config , only : gcr_reduction_target

        use solution , only : nq, q, res, compute_primative_jacobian, dtau, inv_ncells, compute_primative_jacobian

        use grid , only : ncells, cell

        use residual , only : compute_residual

        implicit none

        real(p2), dimension(nq,ncells), intent(in) :: sol_update
        real(p2),                       intent(in) :: gcr_res_rms
        integer,                        intent(out):: iostat

        real(p2), dimension(:,:), pointer   :: q_n
        real(p2), dimension(:,:), pointer   :: r_0
        real(p2)                            :: residual_reduct_target
        real(p2)                            :: Rtau_rms, R0_rms
        real(p2)                            :: delQ_rms, Qn_rms
        real(p2)                            :: delQ_norm
        real(p2)                            :: f_0, f_1, g_1     ! terms in the optimization equation 22
        real(p2), dimension(nq,ncells)      :: frechet_deriv
        real(p2)                            :: ur_opt

        integer :: icell

        q_n => q
        r_0 => res

        nullify(q,res)

        allocate(  q(nq,ncells))
        allocate(res(nq,ncells))

        q = q_n + sol_update

        call compute_residual

        do icell = 1,ncells
            res(:,icell) = res(:,icell) + cell(icell)%vol/dtau(icell) *  &
                            matmul( compute_primative_jacobian(q(:,icell)) , sol_update(:,icell) )
        end do  

        residual_reduct_target = half * (one + gcr_reduction_target)

        Rtau_rms = rms(nq,ncells,res,inv_ncells)
        R0_rms   = rms(nq,ncells,r_0,inv_ncells)

        if ( Rtau_rms / R0_rms < residual_reduct_target ) then
            ! The reduction is acceptable without underrelaxation
            ! Additionally q and res have already been updated.
            deallocate(q_n)
            deallocate(r_0)
            iostat = GCR_SUCCESS
            return
        endif

        ! If the change isn't succesful first we will check if the change is comperable to the computer percision
        delQ_rms = rms(nq,ncells,sol_update,inv_ncells)
        Qn_rms   = rms(nq,ncells,q_n       ,inv_ncells)

        write(*,"(a,es12.6)") "delQ_rms/Qn_rms = ", delQ_rms / Qn_rms 
        if ( delQ_rms / Qn_rms < 1.0e-12_p2) then
            if ( Rtau_rms / R0_rms < one) then
                ! Any reduction will be considered a success at this point, but we don't want to make the CFL any bigger
                deallocate(q_n)
                deallocate(r_0)
                iostat = GCR_CFL_FREEZE
                return
            endif
        endif

        ! The residual did not reduce so now we will apply an underrelaxation factor to minimize the residual
        f_0 = R0_rms
        f_1 = Rtau_rms

        ! We need to compute the frechet derivative again for g_1
        deallocate(q)
        deallocate(res)
        q   => q_n
        res => r_0
        nullify(q_n,r_0)

        ! delQ_norm = l2norm(nq,ncells,sol_update)
        ! call compute_frechet(sol_update,delQ_norm,Qn_rms,frechet_deriv)
        
        ! ! We will temporarily reuse the res vector to save memory space
        ! do icell = 1,ncells
        !     ! EQ 21 from FUN 3D paper where omega = 1
        !     res(:,icell) = res(:,icell) + sol_update(:,icell) * cell(icell)%vol/dtau(icell) + frechet_deriv(:,icell)
        ! end do
        ! g_1 = rms(nq,ncells,res,inv_ncells)
        
        ! We already have the g_1 term from our last projection
        g_1 = gcr_res_rms

        ! Minimize the quadratic function a*ur**2 + b*ur + c
        ! where c = f_0, b = g_1-c, and a = f_1-b-c.
        ur_opt = ( g_1 - f_0 ) / (two * (f_1 - g_1) )
        ! Adding a bound to the ur factor 0 <= ur <= 1
        ur_opt = min(max(ur_opt,zero) , one)

        q = q + ur_opt * sol_update

        ! Check convergence of the updated solution
        call compute_residual

        ! Compute R_tau
        do icell = 1,ncells
            res(:,icell) = res(:,icell) + cell(icell)%vol/dtau(icell) * &
                        matmul( compute_primative_jacobian(q(:,icell)) , sol_update(:,icell) )
        end do  

        Rtau_rms = rms(nq,ncells,res,inv_ncells)

        if ( Rtau_rms / R0_rms < residual_reduct_target .OR. delQ_rms / Qn_rms < 1.0e-12_p2) then
            ! q and res have already been updated and the residual has reduced.
            iostat = GCR_CFL_FREEZE
            return
        endif

        ! If we've made it this far we failed :(
        ! Undo the solution update and report failure
        q = q - ur_opt * sol_update
        iostat = GCR_STALL

    end subroutine gcr_nl_control

    subroutine gcr_CFL_control(gcr_status)

        use common          , only : p2, two
        
        use config          , only : CFL, CFL_max, CFL_min

        use grid            , only : ncells, cell

        use solution        , only : jac, q, compute_primative_jacobian, nq, dtau, compute_local_time_step_dtau, CFL_used

        use direct_solve    , only : gewp_solve

        implicit none

        integer, intent(in) :: gcr_status ! not actually used for now...

        real(p2), dimension(nq,nq) :: preconditioner

        integer :: icell, k, j
        integer :: idestat


        if (gcr_status == GCR_SUCCESS) then
            if (gcr_verbosity >= 3) then
                write(*,*) "Successful iteration:"
                write(*,*) "CFL used:   ", CFL
            endif
            CFL_used = CFL
            CFL = min(CFL * two, CFL_max)
            if (gcr_verbosity >= 3) then
                write(*,*) "CFL update: ", CFL
            endif
        elseif (gcr_status == GCR_CFL_FREEZE) then
            ! Nothing to actually do here
        else ! fail
            do icell = 1,ncells
                preconditioner = compute_primative_jacobian(q(:,icell))

                ! we want to remove the pseudo transient term so that we can add it with a different CFL
                jac(icell)%diag = jac(icell)%diag - (cell(icell)%vol/dtau(icell))*preconditioner
                
            end do

            if (gcr_verbosity >= 3) then
                write(*,*) "Stall detected:"
                write(*,*) "CFL old: ", CFL
            endif
            CFL = max(CFL / 10.0_p2, CFL_min)
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

        endif

    end subroutine gcr_CFL_control

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
        select case(nq)
            ! Loop unrolling
            case(5)
                do i = 1,ncells
                    rms = rms + vector(1,i)**2 &
                              + vector(2,i)**2 &
                              + vector(3,i)**2 &
                              + vector(4,i)**2 &
                              + vector(5,i)**2
                end do
            case default
                do i = 1,ncells
                    do j = 1,nq
                        rms = rms + vector(j,i)**2
                    end do
                end do
        end select

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

        select case(nq)
            ! Loop unrolling
            case(5)
                do i = 1,ncells
                    l2norm = l2norm + vector(1,i)**2 &
                                    + vector(2,i)**2 &
                                    + vector(3,i)**2 &
                                    + vector(4,i)**2 &
                                    + vector(5,i)**2
                end do
            case default
                do i = 1,ncells
                    do j = 1,nq
                        l2norm = l2norm + vector(j,i)**2
                    end do
                end do
        end select

        l2norm = sqrt(l2norm)
    end function

    pure function inner_product(nq,ncells,vector1,vector2)

        ! Function for computing the inner_product of block vectors

        use common , only : p2,zero 

        implicit none

        integer, intent(in)                 :: nq, ncells
        real(p2),dimension(:,:), intent(in) :: vector1, vector2
        real(p2)                            :: inner_product

        integer :: i,j

        inner_product = zero

        select case(nq)
            ! Loop unrolling
            case(5)
                do i = 1,ncells
                    inner_product = inner_product + vector1(1,i) * vector2(1,i) &
                                                  + vector1(2,i) * vector2(2,i) &
                                                  + vector1(3,i) * vector2(3,i) &
                                                  + vector1(4,i) * vector2(4,i) &
                                                  + vector1(5,i) * vector2(5,i) 
                end do
            case default
                do i = 1,ncells
                    do j = 1,nq
                        inner_product = inner_product + vector1(j,i) * vector2(j,i)
                    end do
                end do
        end select

    end function

end module gcr