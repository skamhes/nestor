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

        use config              , only : solver_type, lrelax_sweeps, lrelax_tolerance, smoother

        use grid                , only : ncells, cell

        use solution            , only : jacobian_type

        use sparse_matrix       , only : build_A_BCSM

        ! use gauss_seidel

        use algebraic_multigird , only : UP, DOWN

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
        integer                     :: level = 1
        integer                     :: direction

        call build_A_BCSM(ncells,cell,jacobian_block,V,C,R,nnz=nnz)

        call build_Dinv_array(ncells,num_eq,jacobian_block,Dinv)

        direction = UP 

        call linear_sweeps(num_eq,V,C,R,nnz,residual,Dinv,level,direction,correction,os)

    end subroutine linear_relaxation_block

    subroutine linear_sweeps(num_eq,V,C,R,nnz,res,Dinv,level,direction,correction,stat)

        use common          , only : p2, zero, one, half, two, my_eps

        use config          , only : lrelax_sweeps, solver_type, lrelax_tolerance, smoother, &
                                     use_amg, max_amg_levels

        use solution        , only : lrelax_sweeps_actual, lrelax_roc, roc

        use sparse_matrix   , only : sparseblock_times_vectorblock, build_A_BCSM

        use gauss_seidel    , only : FORWARD, BACKWARD, gauss_seidel_sweep

        use grid            , only : ncells

        implicit none
        
        integer, intent(in)                                 :: num_eq
        real(p2), dimension(:,:,:),allocatable, intent(in)  :: V    ! Values of A
        integer , dimension(:), allocatable, intent(in)     :: C    ! Column index of A
        integer , dimension(ncells+1), intent(in)           :: R    ! Start index of A
        integer , intent(in)                                :: nnz  ! # non-zero terms in A
        real(p2), dimension(num_eq,ncells), intent(in)           :: res  ! RHS (= -b)
        real(p2), dimension(num_eq,num_eq,ncells), intent(in)         :: Dinv ! Inverse of A(i,i)
        integer, intent(in)                                 :: level
        integer, intent(inout)                              :: direction
        
        real(p2), dimension(5,ncells), intent(out)          :: correction
        integer,                       intent(out)          :: stat ! Return 

        ! Residual norms(L1,L2,Linf)
        real(p2), dimension(5)    :: linear_res_norm, linear_res_norm_init

        integer                     :: ii

        ! Under-relaxation parameter
        real(p2) :: omega_lrelax

        ! Initialize some variables
        omega_lrelax         = one
        
        ! Initialize the correction
        correction = zero

        relax : do ii =  1,lrelax_sweeps
            ! Perform Sweep
            if ( trim(smoother) =='gs' ) then
                call gauss_seidel_sweep(num_eq,res,V,R,C,Dinv,nnz,omega_lrelax,correction, linear_res_norm)
            else
                write(*,*) " Sorry, only 'gs' is available at the moment..."
                write(*,*) " Set lrelax_scheme = 'gs', and try again. Stop."
                stop
            endif

            ! Check residual for convergence 
            ! After the first iteration
            if (ii == 1) then
                linear_res_norm_init = linear_res_norm
                ! check for convergence to machine tolerance
                if ( maxval(linear_res_norm(:)) < my_eps ) then
                    lrelax_sweeps_actual = ii
                    stat = RELAX_SUCCESS
                    exit relax
                endif
            else
                roc = maxval(linear_res_norm(1:5)/linear_res_norm_init(1:5))
                if (roc < lrelax_tolerance) then
                    ! if converged
                    lrelax_sweeps_actual = ii
                    stat = RELAX_SUCCESS
                    exit relax
                elseif ( roc > DIVERGENCE_TOLERANCE ) then
                    ! residual has diverged
                    lrelax_sweeps_actual = -1
                    stat = RELAX_FAIL_DIVERGE
                    exit relax
                endif
            endif
        end do relax

        ! ******
        ! 
        ! Algebraic multigrid goes here
        ! 
        ! ******

        if ( roc > lrelax_tolerance .AND. roc < DIVERGENCE_TOLERANCE ) then 
            ! If we make it here the sweeps did not converge
            stat = RELAX_FAIL_STALL
        endif
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
end module linear_solver