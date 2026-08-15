module linear_solver

    use common          , only : p2
    implicit none

    ! This will allow the function to interface with scalar and block matrices
    interface linear_relaxation
        module procedure linear_relaxation_block
        module procedure linear_relaxation_scalar
    end interface linear_relaxation

    interface linear_sweeps
        module procedure linear_sweeps_block
        module procedure linear_sweeps_scalar
    end interface linear_sweeps

    interface multilevel_cycle
        module procedure multilevel_cycle_block
        module procedure multilevel_cycle_scalar
    end interface multilevel_cycle

    public :: RELAX_SUCCESS
    public :: RELAX_FAIL_DIVERGE, RELAX_FAIL_STALL
    integer, parameter  :: RELAX_SUCCESS        = 0
    integer, parameter  :: RELAX_FAIL_DIVERGE   = 1
    integer, parameter  :: RELAX_FAIL_STALL     = 2
    private :: DIVERGENCE_TOLERANCE
    real(p2), parameter :: DIVERGENCE_TOLERANCE = 1e+04

    contains

    ! Solution to the block linear system Ax=b where:
    ! A = jacobian_block is a block matrix with NQxNQ blocks
    ! b = residual block vector with 1xNQ blocks
    ! x = correction is the solution to x=A^(-1)b (also a block vector)
    ! num_eq is the size of the blocks
    subroutine linear_relaxation_block(num_eq, V, Dinv,residual,correction,iostat)

        use common              , only : p2

        use config              , only : solver_type, lrelax_sweeps, lrelax_tolerance, smoother, amg_cycle, line_implicit

        use grid                , only : ncells, cell

        use solution_vars       , only : C, R, nnz

        ! use gauss_seidel

        use algebraic_multigird , only : convert_amg_c_to_i

        implicit none

        integer,                             intent( in) :: num_eq
        real(p2),          dimension(:,:,:), intent( in) :: V
        real(p2),          dimension(:,:,:), intent( in) :: Dinv
        real(p2),            dimension(:,:), intent( in) :: residual
                 
        real(p2),            dimension(:,:), intent(out) :: correction
        integer,                             intent(out) :: iostat

        integer                     :: cycle_type

        cycle_type = convert_amg_c_to_i(amg_cycle)

        if (line_implicit) then
            call li_cycle_block(ncells, num_eq, V,C,R,residual, Dinv, correction, iostat)
        else
            call multilevel_cycle_block(ncells,num_eq,V,C,R,residual,Dinv,cycle_type,.false.,correction,iostat)
        endif

    end subroutine linear_relaxation_block

    subroutine multilevel_cycle_block(ncells,num_eq,V,C,R,res,Dinv,cycle_type,keep_A,correction,stat)

        use common          , only : p2, zero

        use config          , only : lrelax_tolerance, max_amg_cycles

        use solution_vars   , only : lrelax_sweeps_actual, lrelax_roc, roc

        use algebraic_multigird , only : amg_level_block_type, build_amg_struct_block, amg_destroy

        implicit none
        
        integer,                            intent(in)      :: ncells
        integer,                            intent(in)      :: num_eq
        real(p2), dimension(:,:,:), target, intent(in)      :: V    ! Values of A
        integer , dimension(:),     target, intent(in)      :: C    ! Column index of A
        integer , dimension(:),     target, intent(in)      :: R    ! Start index of A
        real(p2), dimension(:,:),           intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(:,:,:), target, intent(in)      :: Dinv ! Inverse of A(i,i)
        integer,                            intent(in)      :: cycle_type
        logical,                            intent(in)      :: keep_A
        
        real(p2), dimension(:,:),           intent(out)     :: correction
        integer,                            intent(out)     :: stat ! Return 

        type(amg_level_block_type) :: base_level

        ! Residual norms(L1,L2,Linf)
        real(p2), dimension(num_eq)    :: linear_res_norm, linear_res_norm_init

        integer                     :: icell, icycle

        ! Compute initial linear residual norm
        linear_res_norm_init = zero
        do icell = 1,ncells
            linear_res_norm_init = linear_res_norm_init + abs(res(:,icell))
        end do

        ! Initialize the correction
        correction = zero

        base_level = build_amg_struct_block(V,C,R,Dinv)
        base_level%ncells = ncells

        amg_cycles : do icycle = 1,max_amg_cycles

            call linear_sweeps(num_eq,base_level,res,cycle_type,1,correction,linear_res_norm,stat)
            
            ! Check for convergence
            roc = maxval(linear_res_norm(1:5)/linear_res_norm_init(1:5))
            if (roc < lrelax_tolerance) then
                ! if converged
                lrelax_sweeps_actual = icycle
                stat = RELAX_SUCCESS
                exit amg_cycles
            elseif ( roc > DIVERGENCE_TOLERANCE ) then
                ! residual has diverged
                lrelax_sweeps_actual = -1
                stat = RELAX_FAIL_DIVERGE
                exit amg_cycles
            endif
        end do amg_cycles


        if ( roc > lrelax_tolerance .AND. roc < DIVERGENCE_TOLERANCE ) then 
            ! If we make it here the sweeps did not converge
            stat = RELAX_FAIL_STALL
        endif
        lrelax_roc = roc

        ! Destroy the amg levels
        if (keep_A) then
            if ( associated(base_level%V) )     nullify(base_level%V)
            if ( associated(base_level%C) )     nullify(base_level%C)
            if ( associated(base_level%R) )     nullify(base_level%R)
            if ( associated(base_level%Dinv) )  nullify(base_level%Dinv)
        endif
        call amg_destroy(base_level)

    end subroutine multilevel_cycle_block

    recursive subroutine linear_sweeps_block(num_eq,solve_level,res,cycle_type,level,correction,l1_res_norm,stat)

        use common          , only : p2, zero, one

        use config          , only : lrelax_sweeps, solver_type, lrelax_tolerance, &
                                     use_amg, max_amg_levels, pre_sweeps, post_sweeps, min_amg_blcoks

        use utils           , only : ismoother, SMOOTH_GS
        
        use sparse_block_matrix   , only : sparseblock_times_vectorblock

        use gauss_seidel    , only : gauss_seidel_sweep

        use algebraic_multigird , only : algebraic_multigrid_prolong, convert_amg_c_to_i, & !,algebraic_multigrid_restrict
                                         AMG_F, AMG_W, AMG_V, amg_restric_rs, amg_level_block_type, build_amg_struct_block

        use ruge_stuben , only : rs_agglom

        implicit none
        
        ! integer, intent(in)                                 :: ncells
        integer,                        intent(in)      :: num_eq
        type(amg_level_block_type), target,   intent(inout)   :: solve_level
        real(p2), dimension(:,:),       intent(in)      :: res  ! RHS (= -b)
        integer,                        intent(in)      :: level
        integer,                        intent(in)      :: cycle_type
        
        real(p2), dimension(:,:),       intent(inout)   :: correction

        integer,                        intent(out)     :: stat ! Return 
        real(p2), dimension(:) ,        intent(out)     :: l1_res_norm

        ! Residual norms(L1,L2,Linf)
        real(p2), dimension(num_eq)    :: linear_res_norm

        integer                     :: isweep

        ! Under-relaxation parameter
        real(p2) :: omega_lrelax

        type(amg_level_block_type), pointer :: coarse_level

        ! Variables to be computed by restricted multigrid
        integer                             :: nnz_restrict     ! number of non-zero entries in restricted matrix
        real(p2), dimension(:,:),   pointer :: restricted_res   ! R*(A*x + b), dimension (nq x ngroups)
        real(p2), dimension(:,:),   pointer :: restricted_correction ! Correction of the restricted linear system 

        ! Ensure the RAP pointers are not in an undefined state:
        nullify(restricted_correction)
        nullify(restricted_res)

        ! Initialize some variables
        omega_lrelax         = one

        pre_sweeps  = max(0,pre_sweeps) ! I use the sweeps to calculate the residual norm
        post_sweeps = max(1,pre_sweeps + post_sweeps) ! we need at least one total sweep

        ! Pre-sweeps
        pre_sweep_loop : do isweep = 1,pre_sweeps
            ! Perform sweeps
            select case(ismoother)
            case(SMOOTH_GS)
                call gauss_seidel_sweep(num_eq,solve_level%ncells,res,solve_level%V,solve_level%C,solve_level%R,solve_level%Dinv, &
                                        omega_lrelax,correction, linear_res_norm)
            case default
                write(*,*) " Sorry, only 'gs' is available at the moment..."
                write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop. linear_solver.f90"
                stop
            end select

        end do pre_sweep_loop

        ! Restrict system
        if (use_amg .and. level < max_amg_levels .and. solve_level%ncells > min_amg_blcoks) then
            ! Restrict the current linear system
            call amg_restric_rs(num_eq,correction,solve_level,res,level, & ! input
                                            nnz_restrict,coarse_level,restricted_res) ! output
            
            ! prepare the restricted correction
            allocate(restricted_correction(num_eq,coarse_level%ncells))

            ! Initialize the restricted correction
            restricted_correction = zero 

            ! Recursively call the linear solver on the restricted system
            call linear_sweeps_block(num_eq,coarse_level,restricted_res,cycle_type,level + 1,& ! input
                            restricted_correction, & ! inout
                            linear_res_norm,stat)    ! Output
            
            ! For the F and W cycles, recursively call the AMG solver again.
            select case(cycle_type)
            case(AMG_F)
                ! Recursively call the linear solver
                call linear_sweeps_block(num_eq,coarse_level,restricted_res,AMG_V,level + 1, & ! input
                                    restricted_correction, & ! inout
                                    linear_res_norm,stat) ! Output
            case(AMG_W)
                ! Recursively call the linear solver
                call linear_sweeps_block(num_eq,coarse_level,restricted_res,cycle_type,level + 1,&!input
                                restricted_correction,& ! inout
                                linear_res_norm,stat) ! Output
            end select
            
            ! Prolong the restricted correction back to its original length
            call algebraic_multigrid_prolong(solve_level%ncells,coarse_level%prolongC,restricted_correction,correction)

        end if

        ! Post-sweeps
        post_sweep_loop : do isweep = 1,post_sweeps
            ! Perform sweeps
            select case(ismoother)
            case(SMOOTH_GS)
                call gauss_seidel_sweep(num_eq,solve_level%ncells,res,solve_level%V,solve_level%C,solve_level%R,solve_level%Dinv, &
                                        omega_lrelax,correction, linear_res_norm)
                case default
                    write(*,*) " Sorry, only 'gs' is available at the moment..."
                    write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop. linear_solver.f90"
                    stop
                end select

        end do post_sweep_loop
        
        ! Make sure all of the allocated arrays are deallocated
        if (associated(restricted_correction))  deallocate(restricted_correction)
        if (associated(restricted_res))         deallocate(restricted_res)

        l1_res_norm = linear_res_norm

    end subroutine linear_sweeps_block

    subroutine li_cycle_block(ncells, num_eq, V,C,R,res, Dinv, correction, stat)

        use common , only : p2, zero

        use config , only : lrelax_sweeps, lrelax_tolerance

        use limplicit , only : lines, nlines

        use solution_vars , only : Rline, iRow, inv_ncells, roc, lrelax_sweeps_actual

        ! use solution_vars   , only : , lrelax_roc


        implicit none

        integer,                            intent(in)      :: ncells
        integer,                            intent(in)      :: num_eq
        real(p2), dimension(:,:,:), target, intent(in)      :: V    ! Values of A
        integer , dimension(:),     target, intent(in)      :: C    ! Column index of A
        integer , dimension(:),     target, intent(in)      :: R    ! Start index of A
        real(p2), dimension(:,:),           intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(:,:,:), target, intent(in)      :: Dinv ! Inverse of A(i,i)
        
        real(p2), dimension(:,:),           intent(out)     :: correction
        integer,                            intent(out)     :: stat ! Return 
        
        integer :: i, ci, k, icell, isweep
        integer :: igs, igse ! start and end pointers for the Gauss-Seidel Sweeps 
        real(p2), dimension(num_eq) :: b, linear_res 
        real(p2), dimension(num_eq) :: linear_res_norm, linear_res_norm_init

        ! Initialize the correction
        correction = zero

        linear_res_norm = zero

        linear_res_norm_init = zero
        do icell = 1,ncells
            linear_res_norm_init = linear_res_norm_init + abs(res(:,icell))
        end do

        sloop : do isweep = 1,lrelax_sweeps
            ! Line sweeps
            do i = 1,nlines
                call thomas_sweep_block(lines(i)%ncells, lines(i)%lcells, num_eq, V, C, R, Rline(2*i-1:2*i+1), Dinv, res, & 
                                        correction, linear_res_norm, stat)
            end do

            ! Gauss-Seidel Sweeps
            igs  = Rline(2*nlines + 1)
            igse = Rline(2*(nlines+1))
            do i = igs,igse-1
                ci = iRow(i)
                ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]
                b = -res(:,ci)
                gs_row_loop : do k = R(i),(R(i+1)-1)
                    ! Add RHS from off diagonal terms and du (du = zero to start and will be updated as we go)
                    if ( C(k) .NE. ci) then
                        b = b - matmul(V(:,:,k),correction(:,C(k)))
                    end if
                end do gs_row_loop
                ! ! Update du by the GS relaxation:
                !
                ! e.g., for 3 nghbrs, perform the relaxation in the form:
                !
                !                     diagonal block        sum of off-diagonal block contributions
                !       dUj = omega*{ [V/dtj+dR/dUj]^{-1}*(-[dRj/dU1]*dU1 -[dRj/dU2]*dU2 -[dRj/dU3]*dU3 -Res_j) - dUj }
                linear_res = matmul(Dinv(:,:,ci), b) - correction(:,ci)
                correction(:,ci) = correction(:,ci) +  linear_res
                linear_res_norm(:) = linear_res_norm(:) + abs(linear_res)

            end do

            linear_res_norm = linear_res_norm / real(ncells, p2)

            ! Check for convergence
            roc = maxval(linear_res_norm(1:5)/linear_res_norm_init(1:5))
            if (roc < lrelax_tolerance) then
                ! if converged
                lrelax_sweeps_actual = isweep
                stat = RELAX_SUCCESS
                exit sloop
            elseif ( roc > DIVERGENCE_TOLERANCE ) then
                ! residual has diverged
                lrelax_sweeps_actual = -1
                stat = RELAX_FAIL_DIVERGE
                exit sloop
            endif
        end do sloop

    end subroutine li_cycle_block

    subroutine thomas_sweep_block(nc, lcells, neq, V, C, R, Rline, Dinv, res, correction, linear_res, stat)

        use common , only : p2

        use direct_solve        , only : gewp_solve

        implicit none
        
        integer,                            intent(in)      :: nc     ! number of cells in the given line
        integer,  dimension(:),             intent(in)      :: lcells ! Array of cells in the line
        integer,                            intent(in)      :: neq
        real(p2), dimension(:,:,:), target, intent(in)      :: V      ! Values of A
        integer , dimension(:),     target, intent(in)      :: C      ! Column index of A
        integer , dimension(:),     target, intent(in)      :: R      ! Start index of each row in A
        integer , dimension(3),     target, intent(in)      :: Rline  ! Start index of R for each line block in A
        real(p2), dimension(:,:),           intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(:,:,:), target, intent(in)      :: Dinv ! Inverse of A(i,i)
        
        real(p2), dimension(:,:),           intent(inout)   :: correction
        real(p2), dimension(neq),           intent(inout)   :: linear_res
        integer,                            intent(out)     :: stat ! Return 

        real(p2), dimension(neq,nc) :: rhs ! scratch space. We can't clobber res or correction (at first)
        real(p2), dimension(neq,neq,nc) :: deltai

        integer :: po, pl, pn ! pointers to the off line block, line block, and the next block
        
        integer :: i, j, k, jj
        integer :: ci, cj

        real(p2), dimension(neq,neq) :: l, d, u, um1, di ! 3 tridiagonal blocks
        real(p2), dimension(neq)     :: new_corr, lres
        
        po = Rline(1)
        pl = Rline(2)
        pn = Rline(3)

        ! First add the off loop blocks
        k = 0
        do i = po,pl-1
            k = k+1
            ci = lcells(k)
            rhs(:,k) = - res(:,ci)
            do j = R(i),R(i+1)-1
                cj = C(j)
                rhs(:,k) = rhs(:,k) - matmul(V(:,:,j),correction(:,cj))
            end do
        end do

        ! Now perform the Thomas Algorithm
        ! This implements the algorithm in https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Method
        ! but implemented with blocks
        
        ! The first and last rows have special treatment:
        deltai(:,:,1) = matmul(Dinv(:,:,lcells(1)), V(:,:,R(pl)+1))
        rhs(:,1)      = matmul(Dinv(:,:,lcells(1)), rhs(:,1))
        
        j = 1 ! local cell counter
        do i = pl+1,pn-2
            j = j+1
            l   = V(:,:,R(i)  )
            d   = V(:,:,R(i)+1)
            u   = V(:,:,R(i)+2)

            di = d - matmul(l,deltai(:,:,j-1)) !MM
            call gewp_solve( di, neq  , deltai(:,:,j), stat) ! fortran doesn't let you alias variables
            !  Report errors
            if (stat/=0) then
                write(*,*) " Error in inverting the diagonal block... Stop"
                write(*,*) "  Cell number = ", lcells(j)
                do k = 1, neq
                    write(*,'(12(es8.1))') ( di(k,jj), jj=1,5 )
                end do
                stop
            endif

            rhs(:,j)      = rhs(:,j) - matmul(l,rhs(:,j-1)) !MV
            rhs(:,j)      = matmul(deltai(:,:,j),rhs(:,j))  !MV
            deltai(:,:,j) = matmul(deltai(:,:,j),u)         !MM
        end do

        j = j+1
        l   = V(:,:,R(pn-1)  )
        d   = V(:,:,R(pn-1)+1)
        di = d - matmul(l,deltai(:,:,j-1))
        call gewp_solve( di, neq  , deltai(:,:,j), stat) ! fortran doesn't let you alias variables
        rhs(:,j)      = rhs(:,j) - matmul(l,rhs(:,j-1))
        rhs(:,j)      = matmul(deltai(:,:,j),rhs(:,j))


        ! now we back substitute to update the correction
        new_corr = rhs(:,j)
        linear_res = linear_res + abs(new_corr - correction(:,lcells(j)))
        correction(:,lcells(j)) = new_corr
        do i = pn-1,pl+1,-1 ! loop backwards
            j = j - 1
            u = V(:,:,C(R(i-1)))
            new_corr = rhs(:,j) - matmul(deltai(:,:,j),correction(:,lcells(j+1))) !MV
            linear_res = linear_res + abs(new_corr - correction(:,lcells(j)))
            correction(:,lcells(j)) = new_corr
        end do


    end subroutine thomas_sweep_block

    subroutine linear_relaxation_scalar(V,Dinv,residual,correction,iostat)

        use common              , only : p2

        use config              , only : solver_type, lrelax_sweeps, lrelax_tolerance, smoother, amg_cycle, line_implicit

        use grid                , only : ncells, cell

        use solution_vars       , only : C, R, nnz

        use algebraic_multigird , only : convert_amg_c_to_i

        implicit none

        real(p2), dimension(:), intent( in) :: V
        real(p2), dimension(:), intent( in) :: Dinv
        real(p2), dimension(:), intent( in) :: residual
                 
        real(p2), dimension(:), intent(out) :: correction
        integer,                intent(out) :: iostat

        integer                             :: cycle_type

        cycle_type = convert_amg_c_to_i(amg_cycle)
        if (line_implicit) then
            call li_cycle_scalar(ncells, V,C,R,residual, Dinv, correction, iostat)
        else
            call multilevel_cycle(ncells,V,C,R,residual,Dinv,cycle_type,.false.,correction,iostat)
        endif
    end subroutine linear_relaxation_scalar

    subroutine multilevel_cycle_scalar(ncells,V,C,R,res,Dinv,cycle_type,keep_A,correction,stat)

        use common          , only : p2, zero

        use config          , only : lrelax_tolerance, max_amg_cycles

        use solution_vars   , only : lrelax_sweeps_actual, lrelax_roc, roc

        use algebraic_multigird , only : amg_level_scalar_type, build_amg_struct_scalar, amg_destroy

        implicit none
        
        integer,                        intent(in)      :: ncells
        real(p2), dimension(:), target, intent(in)      :: V    ! Values of A
        integer , dimension(:), target, intent(in)      :: C    ! Column index of A
        integer , dimension(:), target, intent(in)      :: R    ! Start index of A
        real(p2), dimension(:),         intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(:), target, intent(in)      :: Dinv ! Inverse of A(i,i)
        integer,                        intent(in)      :: cycle_type
        logical,                        intent(in)      :: keep_A
        
        real(p2), dimension(:),         intent(out)     :: correction
        integer,                        intent(out)     :: stat ! Return 

        type(amg_level_scalar_type) :: base_level

        ! Residual norms(L1,L2,Linf)
        real(p2)    :: linear_res_norm, linear_res_norm_init

        integer     :: icell, icycle

        ! Compute initial linear residual norm
        linear_res_norm_init = zero
        do icell = 1,ncells
            linear_res_norm_init = linear_res_norm_init + abs(res(icell))
        end do

        ! Initialize the correction
        correction = zero

        base_level = build_amg_struct_scalar(V,C,R,Dinv)
        base_level%ncells = ncells

        amg_cycles : do icycle = 1,max_amg_cycles

            call linear_sweeps_scalar(base_level,res,cycle_type,1,correction,linear_res_norm,stat)
            
            ! Check for convergence
            roc = linear_res_norm/linear_res_norm_init
            if (roc < lrelax_tolerance) then
                ! if converged
                lrelax_sweeps_actual = icycle
                stat = RELAX_SUCCESS
                exit amg_cycles
            elseif ( roc > DIVERGENCE_TOLERANCE ) then
                ! residual has diverged
                lrelax_sweeps_actual = -1
                stat = RELAX_FAIL_DIVERGE
                exit amg_cycles
            endif
        end do amg_cycles


        if ( roc > lrelax_tolerance .AND. roc < DIVERGENCE_TOLERANCE ) then 
            ! If we make it here the sweeps did not converge
            stat = RELAX_FAIL_STALL
        endif
        lrelax_roc = roc

        ! Destroy the amg levels
        if (keep_A) then
            if ( associated(base_level%V) )     nullify(base_level%V)
            if ( associated(base_level%C) )     nullify(base_level%C)
            if ( associated(base_level%R) )     nullify(base_level%R)
            if ( associated(base_level%Dinv) )  nullify(base_level%Dinv)
        endif
        call amg_destroy(base_level)

    end subroutine multilevel_cycle_scalar

    recursive subroutine linear_sweeps_scalar(solve_level,res,cycle_type,level,correction,l1_res_norm,stat)

        use common          , only : p2, zero, one

        use config          , only : lrelax_sweeps, solver_type, lrelax_tolerance, &
                                     use_amg, max_amg_levels, pre_sweeps, post_sweeps, min_amg_blcoks

        use utils           , only : ismoother, SMOOTH_GS
        
        use sparse_scalar_matrix   , only : sparseMat_times_vector

        use gauss_seidel    , only : gauss_seidel_sweep

        use algebraic_multigird , only : algebraic_multigrid_prolong, convert_amg_c_to_i, &
                                         AMG_F, AMG_W, AMG_V, amg_restric_rs, amg_level_scalar_type, build_amg_struct_block

        use ruge_stuben , only : rs_agglom

        implicit none
        
        type(amg_level_scalar_type), target, intent(inout)   :: solve_level
        real(p2), dimension(:),              intent(in)      :: res  ! RHS (= -b)
        integer,                             intent(in)      :: level
        integer,                             intent(in)      :: cycle_type
        
        real(p2), dimension(:),              intent(inout)   :: correction

        integer,                             intent(out)     :: stat ! Return 
        real(p2),                            intent(out)     :: l1_res_norm

        ! Residual norms(L1,L2,Linf)
        real(p2)    :: linear_res_norm

        integer                     :: isweep

        ! Under-relaxation parameter
        real(p2) :: omega_lrelax

        type(amg_level_scalar_type), pointer :: coarse_level

        ! Variables to be computed by restricted multigrid
        integer                           :: nnz_restrict     ! number of non-zero entries in restricted matrix
        real(p2), dimension(:),   pointer :: restricted_res   ! R*(A*x + b), dimension (nq x ngroups)
        real(p2), dimension(:),   pointer :: restricted_correction ! Correction of the restricted linear system 

        ! Ensure the RAP pointers are not in an undefined state:
        nullify(restricted_correction)
        nullify(restricted_res)

        ! Initialize some variables
        omega_lrelax         = one

        pre_sweeps  = max(0,pre_sweeps) ! I use the sweeps to calculate the residual norm
        post_sweeps = max(1,pre_sweeps + post_sweeps) ! we need at least one total sweep

        ! Pre-sweeps
        pre_sweep_loop : do isweep = 1,pre_sweeps
            ! Perform sweeps
            select case(ismoother)
            case(SMOOTH_GS)
                call gauss_seidel_sweep(solve_level%ncells,res,solve_level%V,solve_level%C,solve_level%R,solve_level%Dinv, &
                                        omega_lrelax,correction, linear_res_norm)
            case default
                write(*,*) " Sorry, only 'gs' is available at the moment..."
                write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop. linear_solver.f90"
                stop
            end select

        end do pre_sweep_loop

        ! Restrict system
        if (use_amg .and. level < max_amg_levels .and. solve_level%ncells > min_amg_blcoks) then
            ! Restrict the current linear system
            call amg_restric_rs(correction,solve_level,res,level, & ! input
                                            nnz_restrict,coarse_level,restricted_res) ! output
            
            ! prepare the restricted correction
            allocate(restricted_correction(coarse_level%ncells))

            ! Initialize the restricted correction
            restricted_correction = zero 

            ! Recursively call the linear solver on the restricted system
            call linear_sweeps_scalar(coarse_level,restricted_res,cycle_type,level + 1,& ! input
                            restricted_correction, & ! inout
                            linear_res_norm,stat)    ! Output
            
            ! For the F and W cycles, recursively call the AMG solver again.
            select case(cycle_type)
            case(AMG_F)
                ! Recursively call the linear solver
                call linear_sweeps_scalar(coarse_level,restricted_res,AMG_V,level + 1, & ! input
                                    restricted_correction, & ! inout
                                    linear_res_norm,stat) ! Output
            case(AMG_W)
                ! Recursively call the linear solver
                call linear_sweeps_scalar(coarse_level,restricted_res,cycle_type,level + 1,&!input
                                restricted_correction,& ! inout
                                linear_res_norm,stat) ! Output
            end select
            
            ! Prolong the restricted correction back to its original length
            call algebraic_multigrid_prolong(solve_level%ncells,coarse_level%prolongC,restricted_correction,correction)

        end if

        ! Post-sweeps
        post_sweep_loop : do isweep = 1,post_sweeps
            ! Perform sweeps
            select case(ismoother)
            case(SMOOTH_GS)
                call gauss_seidel_sweep(solve_level%ncells,res,solve_level%V,solve_level%C,solve_level%R,solve_level%Dinv, &
                                        omega_lrelax,correction, linear_res_norm)
                case default
                    write(*,*) " Sorry, only 'gs' is available at the moment..."
                    write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop. linear_solver.f90"
                    stop
                end select

        end do post_sweep_loop
        
        ! Make sure all of the allocated arrays are deallocated
        if (associated(restricted_correction))  deallocate(restricted_correction)
        if (associated(restricted_res))         deallocate(restricted_res)

        l1_res_norm = linear_res_norm

    end subroutine linear_sweeps_scalar

    subroutine li_cycle_scalar(ncells, V,C,R,res, Dinv, correction, stat)

        use common , only : p2, zero

        use config , only : lrelax_sweeps, lrelax_tolerance

        use limplicit , only : lines, nlines

        use solution_vars , only : Rline, iRow, inv_ncells, roc, lrelax_sweeps_actual

        ! use solution_vars   , only : , lrelax_roc

        implicit none

        integer,                        intent(in)      :: ncells
        real(p2), dimension(:), target, intent(in)      :: V    ! Values of A
        integer , dimension(:), target, intent(in)      :: C    ! Column index of A
        integer , dimension(:), target, intent(in)      :: R    ! Start index of A
        real(p2), dimension(:),         intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(:), target, intent(in)      :: Dinv ! Inverse of A(i,i)
        
        real(p2), dimension(:),         intent(out)     :: correction
        integer,                        intent(out)     :: stat ! Return 
        
        integer :: i, ci, k, icell, isweep
        integer :: igs, igse ! start and end pointers for the Gauss-Seidel Sweeps 
        real(p2) :: b, linear_res 
        real(p2) :: linear_res_norm, linear_res_norm_init

        ! Initialize the correction
        correction = zero

        linear_res_norm = zero

        linear_res_norm_init = zero
        do icell = 1,ncells
            linear_res_norm_init = linear_res_norm_init + abs(res(icell))
        end do

        sloop : do isweep = 1,lrelax_sweeps
            ! Line sweeps
            do i = 1,nlines
                call thomas_sweep_scalar(lines(i)%ncells, lines(i)%lcells, V, C, R, Rline(2*i-1:2*i+1), Dinv, res, correction, &
                                linear_res_norm, stat)
            end do

            ! Gauss-Seidel Sweeps
            igs  = Rline(2*nlines + 1)
            igse = Rline(2*(nlines+1))
            do i = igs,igse-1
                ci = iRow(i)
                ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]
                b = -res(ci)
                gs_row_loop : do k = R(i),(R(i+1)-1)
                    ! Add RHS from off diagonal terms and du (du = zero to start and will be updated as we go)
                    if ( C(k) .NE. ci) then
                        b = b - V(k)*correction(C(k))
                    end if
                end do gs_row_loop
                ! ! Update du by the GS relaxation:
                !
                ! e.g., for 3 nghbrs, perform the relaxation in the form:
                !
                !                     diagonal block        sum of off-diagonal block contributions
                !       dUj = omega*{ [V/dtj+dR/dUj]^{-1}*(-[dRj/dU1]*dU1 -[dRj/dU2]*dU2 -[dRj/dU3]*dU3 -Res_j) - dUj }
                linear_res = Dinv(ci)*b - correction(ci)
                correction(ci) = correction(ci) +  linear_res
                linear_res_norm = linear_res_norm + abs(linear_res)

            end do

            linear_res_norm = linear_res_norm / real(ncells, p2)

            ! Check for convergence
            roc = linear_res_norm/linear_res_norm_init
            if (roc < lrelax_tolerance) then
                ! if converged
                lrelax_sweeps_actual = isweep
                stat = RELAX_SUCCESS
                exit sloop
            elseif ( roc > DIVERGENCE_TOLERANCE ) then
                ! residual has diverged
                lrelax_sweeps_actual = -1
                stat = RELAX_FAIL_DIVERGE
                exit sloop
            endif
        end do sloop

    end subroutine li_cycle_scalar

    subroutine thomas_sweep_scalar(nc, lcells, V, C, R, Rline, Dinv, res, correction, linear_res, stat)

        use common , only : p2

        use direct_solve        , only : gewp_solve

        implicit none
        
        integer,                intent(in)      :: nc     ! number of cells in the given line
        integer,  dimension(:), intent(in)      :: lcells ! Array of cells in the line
        real(p2), dimension(:), intent(in)      :: V      ! Values of A
        integer , dimension(:), intent(in)      :: C      ! Column index of A
        integer , dimension(:), intent(in)      :: R      ! Start index of each row in A
        integer , dimension(3), intent(in)      :: Rline  ! Start index of R for each line block in A
        real(p2), dimension(:), intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(:), intent(in)      :: Dinv ! Inverse of A(i,i)
        
        real(p2), dimension(:), intent(inout)   :: correction
        real(p2),               intent(inout)   :: linear_res
        integer,                intent(out)     :: stat ! Return 

        real(p2), dimension(nc) :: rhs ! scratch space. We can't clobber res or correction (at first)
        real(p2), dimension(nc) :: deltai

        integer :: po, pl, pn ! pointers to the off line block, line block, and the next block
        
        integer :: i, j, k, jj
        integer :: ci, cj

        real(p2) :: l, d, u, um1, di ! 3 tridiagonal blocks
        real(p2) :: new_corr, lres
        
        po = Rline(1)
        pl = Rline(2)
        pn = Rline(3)

        ! First add the off loop blocks
        k = 0
        do i = po,pl-1
            k = k+1
            ci = lcells(k)
            rhs(k) = - res(ci)
            do j = R(i),R(i+1)-1
                cj = C(j)
                rhs(k) = rhs(k) - V(j)*correction(cj)
            end do
        end do

        ! Now perform the Thomas Algorithm
        ! This implements the algorithm in https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Method
        ! but implemented with blocks
        
        ! The first and last rows have special treatment:
        deltai(1) = Dinv(lcells(1))*V(R(pl)+1)
        rhs(1)    = Dinv(lcells(1))*rhs(1)
        
        j = 1 ! local cell counter
        do i = pl+1,pn-2
            j = j+1
            l   = V(R(i)  )
            d   = V(R(i)+1)
            u   = V(R(i)+2)

            di     = d - l*deltai(j-1)
            deltai(j) = 1 / max(di,1.e-10_p2)
            rhs(j)    = rhs(j) - l* rhs(j-1)
            rhs(j)    = deltai(j) * rhs(j)
            deltai(j) = deltai(j) * u
        end do

        j = j+1
        l  = V(R(pn-1)  )
        d  = V(R(pn-1)+1)
        di = d - l*deltai(j-1)
        deltai(j) = 1 / max(di,1.e-10_p2)
        rhs(j)    = rhs(j) - l*rhs(j-1)
        rhs(j)    = deltai(j) *rhs(j)


        ! now we back substitute to update the correction
        new_corr = rhs(j)
        linear_res = linear_res + abs(new_corr - correction(lcells(j)))
        correction(lcells(j)) = new_corr
        do i = pn-1,pl+1,-1 ! loop backwards
            j = j - 1
            u = V(C(R(i-1)))
            new_corr = rhs(j) - deltai(j)*correction(lcells(j+1))
            linear_res = linear_res + abs(new_corr - correction(lcells(j)))
            correction(lcells(j)) = new_corr
        end do


    end subroutine thomas_sweep_scalar
end module linear_solver
