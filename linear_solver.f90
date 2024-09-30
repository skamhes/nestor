module linear_solver

    use common          , only : p2
    implicit none

    ! This will allow the function to interface with scalar and block matrices
    interface linear_relaxation
        module procedure linear_relaxation_block
    end interface linear_relaxation

    public :: RELAX_SUCCESS
    public :: RELAX_FAIL_DIVERGE, RELAX_FAIL_STALL
    integer, parameter :: RELAX_SUCCESS      = 0
    integer, parameter :: RELAX_FAIL_DIVERGE = 1
    integer, parameter :: RELAX_FAIL_STALL   = 2
    private :: DIVERGENCE_TOLERANCE
    real(p2) :: DIVERGENCE_TOLERANCE = 1e+04

    contains

    ! Solution to the block linear system Ax=b where:
    ! A = jacobian_block is a block matrix with NQxNQ blocks
    ! b = residual block vector with 1xNQ blocks
    ! x = correction is the solution to x=A^(-1)b (also a block vector)
    ! num_eq is the size of the blocks
    subroutine linear_relaxation_block(num_eq,jacobian_block,residual,correction)

        use common              , only : p2, one, &
                                         my_eps

        use config              , only : solver_type, lrelax_sweeps, lrelax_tolerance, smoother, amg_cycle

        use grid                , only : ncells, cell

        use solution            , only : jacobian_type

        ! use gauss_seidel

        use algebraic_multigird , only : UP, DOWN, convert_amg_c_to_i

        implicit none

        integer,                                intent( in) :: num_eq
        type(jacobian_type), dimension(ncells), intent( in) :: jacobian_block
        real(p2), dimension(   num_eq, ncells), intent( in) :: residual

        real(p2), dimension(   num_eq, ncells), intent(out) :: correction

        real(p2), dimension(:,:,:), allocatable :: V   ! Values (5x5 block matrix) plus corresponding index
        integer , dimension(:),     allocatable :: C   ! Column index of each value
        integer , dimension(ncells+1)           :: R   ! Start index of each new row
        integer                                 :: nnz
        real(p2), dimension(5,5,ncells)             :: Dinv

        integer                     :: os
        integer                     :: level
        integer                     :: cycle_type

        call build_A_BCSM(ncells,cell,jacobian_block,V,C,R,nnz=nnz)

        call build_Dinv_array(ncells,num_eq,jacobian_block,Dinv)

        cycle_type = convert_amg_c_to_i(amg_cycle)

        level = 1

        call multilevel_cycle(ncells,num_eq,nnz,V,C,R,residual,Dinv,cycle_type,level,correction,os)

    end subroutine linear_relaxation_block

    subroutine multilevel_cycle(ncells,num_eq,nnz,V,C,R,res,Dinv,cycle_type,level,correction,stat)

        use common          , only : p2, zero

        use config          , only : lrelax_tolerance, max_amg_cycles

        use solution        , only : lrelax_sweeps_actual, lrelax_roc, roc

        implicit none
        
        integer, intent(in)                                 :: ncells
        integer, intent(in)                                 :: num_eq
        integer, intent(in)                                 :: nnz
        real(p2), dimension(num_eq,num_eq,nnz), intent(in)  :: V    ! Values of A
        integer , dimension(nnz)           , intent(in)     :: C    ! Column index of A
        integer , dimension(ncells+1), intent(in)           :: R    ! Start index of A
        real(p2), dimension(num_eq,ncells), intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(num_eq,num_eq,ncells),intent(in):: Dinv ! Inverse of A(i,i)
        integer, intent(inout)                              :: level
        integer, intent(in)                                 :: cycle_type
        
        real(p2), dimension(num_eq,ncells), intent(out)          :: correction
        integer,                       intent(out)          :: stat ! Return 

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

        amg_cycles : do icycle = 1,max_amg_cycles

            call linear_sweeps(ncells,num_eq,nnz,V,C,R,res,Dinv,cycle_type,level,correction,linear_res_norm,stat)
            
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

    end subroutine multilevel_cycle

    recursive subroutine linear_sweeps(ncells,num_eq,nnz,V,C,R,res,Dinv,cycle_type,level,correction,l1_res_norm,stat)

        use common          , only : p2, zero, one, half, two, my_eps

        use config          , only : lrelax_sweeps, solver_type, lrelax_tolerance, smoother, &
                                     use_amg, max_amg_levels, pre_sweeps, post_sweeps, min_amg_blcoks

        use sparse_matrix   , only : sparseblock_times_vectorblock

        use gauss_seidel    , only : FORWARD, BACKWARD, gauss_seidel_sweep

        use algebraic_multigird , only : algebraic_multigrid_prolong, algebraic_multigrid_restrict, UP, DOWN, convert_amg_c_to_i, &
                                         AMG_F, AMG_W, AMG_V

        implicit none
        
        integer, intent(in)                                 :: ncells
        integer, intent(in)                                 :: num_eq
        integer, intent(in)                                 :: nnz
        real(p2), dimension(num_eq,num_eq,nnz), intent(in)  :: V    ! Values of A
        integer , dimension(nnz)           , intent(in)     :: C    ! Column index of A
        integer , dimension(ncells+1), intent(in)           :: R    ! Start index of A
        real(p2), dimension(num_eq,ncells), intent(in)      :: res  ! RHS (= -b)
        real(p2), dimension(num_eq,num_eq,ncells),intent(in):: Dinv ! Inverse of A(i,i)
        integer, intent(in)                                 :: level
        integer, intent(in)                                 :: cycle_type
        
        real(p2), dimension(num_eq,ncells), intent(inout)   :: correction

        integer,                       intent(out)          :: stat ! Return 
        real(p2), dimension(num_eq) , intent(out)           :: l1_res_norm

        ! Residual norms(L1,L2,Linf)
        real(p2), dimension(num_eq)    :: linear_res_norm, linear_res_norm_init

        integer                     :: icell, isweep

        ! Under-relaxation parameter
        real(p2) :: omega_lrelax

        ! Variables to be computed by restricted multigrid
        integer                             :: ngroup           ! Number of groups in level n + 1
        integer                             :: nnz_restrict     ! number of non-zero entries in restricted matrix
        integer,  dimension(ncells)         :: prolongC         ! Column ind. of prolongation matrix. Needed to prolong correction
        real(p2), dimension(:,:,:), pointer :: RAP_V            ! Values of restricted matrix RAP
        integer,  dimension(:),     pointer :: RAP_C            ! Column indices of RAP
        integer,  dimension(:),     pointer :: RAP_R            ! Row indices of RAP (length = ngroups + 1)
        real(p2), dimension(:,:,:), pointer :: RAP_Dinv         ! Inverse of RAP diagonal blocks (dim = nq x nq x ngroups)
        real(p2), dimension(:,:),   pointer :: restricted_res   ! R*(A*x + b), dimension (nq x ngroups)
        real(p2), dimension(:,:),   pointer :: restricted_correction ! Correction of the restricted linear system 

        ! Ensure the RAP pointers are not in an undefined state:
        nullify(RAP_V)
        nullify(RAP_C)
        nullify(RAP_R)
        nullify(RAP_Dinv)
        nullify(restricted_correction)
        nullify(restricted_res)

        ! Initialize some variables
        omega_lrelax         = one

        ! write(*,*) "level: ", level
        
        pre_sweeps  = max(0,pre_sweeps) ! I use the sweeps to calculate the residual norm
        post_sweeps = max(1,pre_sweeps + post_sweeps) ! we need at least one total sweep

        ! Pre-sweeps
        pre_sweep_loop : do isweep = 1,pre_sweeps
            ! Perform sweeps
            if ( trim(smoother) =='gs' ) then
                call gauss_seidel_sweep(num_eq,ncells,res,V,C,R,Dinv,omega_lrelax,correction, linear_res_norm)
            else
                write(*,*) " Sorry, only 'gs' is available at the moment..."
                write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop."
                stop
            endif

        end do pre_sweep_loop

        ! Restrict system
        if (use_amg .and. level < max_amg_levels .and. ncells > min_amg_blcoks) then
            ! Restrict the current linear system
            call algebraic_multigrid_restrict(ncells,num_eq,correction,V,C,R,nnz,res,level, & ! input
                                            ngroup,nnz_restrict,prolongC,RAP_V,RAP_C,RAP_R,RAP_Dinv,restricted_res) ! output
            
            ! prepare the restricted correction
            allocate(restricted_correction(num_eq,ngroup))

            ! Initialize the restricted correction
            restricted_correction = zero 

            ! Recursively call the linear solver
            call linear_sweeps(ngroup,num_eq,nnz_restrict,RAP_V,RAP_C,RAP_R,restricted_res,RAP_Dinv,cycle_type,level + 1,& ! input
                            restricted_correction, & ! inout
                            linear_res_norm,stat)    ! Output
            
            ! For the F and V cycles, recursively call the AMG solver again.
            select case(cycle_type)
            case(AMG_F)
                    ! Recursively call the linear solver
                call linear_sweeps(ngroup,num_eq,nnz_restrict,RAP_V,RAP_C,RAP_R,restricted_res,RAP_Dinv,AMG_V,level + 1, & ! input
                                    restricted_correction, & ! inout
                                    linear_res_norm,stat) ! Output
            case(AMG_W)
                    ! Recursively call the linear solver
                call linear_sweeps(ngroup,num_eq,nnz_restrict,RAP_V,RAP_C,RAP_R,restricted_res,RAP_Dinv,cycle_type,level + 1,&!input
                                restricted_correction,& ! inout
                                linear_res_norm,stat) ! Output
            end select
            
            ! Prolong the restricted correction back to its original length
            call algebraic_multigrid_prolong(ncells,prolongC,restricted_correction,correction)

        end if

        ! Post-sweeps
        post_sweep_loop : do isweep = 1,post_sweeps
            ! Perform sweeps
            if ( trim(smoother) =='gs' ) then
                call gauss_seidel_sweep(num_eq,ncells,res,V,C,R,Dinv,omega_lrelax,correction, linear_res_norm)
            else
                write(*,*) " Sorry, only 'gs' is available at the moment..."
                write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop."
                stop
            endif

        end do post_sweep_loop
        

        ! Make sure all of the allocated arrays are deallocated
        if (associated(RAP_V))                  deallocate(RAP_V)
        if (associated(RAP_C))                  deallocate(RAP_C)
        if (associated(RAP_R))                  deallocate(RAP_R)
        if (associated(RAP_Dinv))               deallocate(RAP_Dinv)
        if (associated(restricted_correction))  deallocate(restricted_correction)
        if (associated(restricted_res))         deallocate(restricted_res)

        l1_res_norm = linear_res_norm

        ! write(*,*) "level: ", level
        ! write(*,*)
    end subroutine linear_sweeps

    subroutine build_Dinv_array(ncells,nq,jac,D_inv)
        ! This subroutine stores just the diagonal blocks of the Dinv matrix since the rest are empty.  As a result, C=R=index
        ! so the other two indices do not need to be stored.
        use common      , only : p2, zero

        use solution    , only : jacobian_type

        implicit none
        integer, intent(in) :: ncells
        integer, intent(in) :: nq
        type(jacobian_type), dimension(ncells), intent(in) :: jac
        
        real(p2), dimension(nq,nq,ncells), INTENT(OUT) :: D_inv
        

        integer :: i
        
        do i = 1,ncells
            D_inv(:,:,i) = jac(i)%diag_inv(:,:)
        end do
    end subroutine build_Dinv_array

    subroutine build_A_BCSM(ncells,cell,jac,V,C,R,nnz)
        ! Takes in a cell structure and a jacobian structure and creates a corresponding A matrix using the Yale meethod:
        ! https://en.wikipedia.org/wiki/Sparse_matrix


        use common      , only : p2

        use grid        , only : cc_data_type

        use solution    , only : jacobian_type

        use sparse_common, only: insertion_sort_index

        implicit none 

        integer, intent(in) :: ncells
        type(cc_data_type), dimension(ncells), intent(in) :: cell
        type(jacobian_type), dimension(ncells), intent(in) :: jac

        real(p2), dimension(:,:,:), allocatable, intent(out)  :: V   ! Values (5x5 block matrix) plus corresponding index
        integer, dimension(:),  allocatable, intent(out)       :: C   ! Column index of each value
        integer,dimension(ncells+1), intent(out) :: R   ! Start index of each new row
        integer, INTENT(OUT), optional :: nnz

        integer :: i, j, length

        R(1) = 1 ! Row 1 starts at 1
        do i = 2,ncells + 1
            R(i) = R(i-1) + 1 + cell(i-1)%nnghbrs ! Start of row(i) = row(i-1) start point + 1 (diagonal term) + # of neighbors
        end do
        nnz = R(ncells+1) - 1 ! number of nonzero cells

        allocate(V(5,5,nnz))
        allocate(C(    nnz))

        do i = 1,ncells
            ! sort the index of the cell neighbors and i and stores them in C:
            call insertion_sort_index( (/ cell(i)%nghbr, i /) , C(R(i) : (R(i+1)-1)) ) 
            length = R(i+1)-R(i)
            do j = R(i),(R(i+1)-1)
                if (length == C(j)) then
                    V(:,:,j) = jac(i)%diag(:,:)
                    C(j) = i
                else
                    V(:,:,j) = jac(i)%off_diag(:,:,C(j))
                    C(j) = cell(i)%nghbr(C(j))
                end if
            end do
        end do
    end subroutine build_A_BCSM

end module linear_solver