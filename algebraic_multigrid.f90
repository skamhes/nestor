module algebraic_multigird

    implicit none

    public :: UP, DOWN
    integer, parameter :: UP    = 1 ! Coarsen
    integer, parameter :: DOWN  = 2 ! Refine

    public :: AMG_F, AMG_V, AMG_W
    integer, parameter :: AMG_F = 1
    integer, parameter :: AMG_W = 2
    integer, parameter :: AMG_V = 3

    contains

    subroutine algebraic_multigrid_restrict_legacy(ncells,nq,phi,V,C,R,nnz,res,level, &
                                                      ngroup,nnz_restrict,prolongC,RAP_V,RAP_C,RAP_R,RAP_Dinv,defect_res)
    
        ! Implementation of algebraic multi grid using additive correction.
        use common                  , only : p2, zero

        use grid                    , only : cc_data_type
        
        use direct_solve            , only : gewp_solve

        use sparse_common           , only : insertion_sort_int, insertion_sort_real


        !use stdlib_sorting, only : sort
        
        implicit none
        ! INPUT
        integer,                                intent(in)  :: ncells ! number of cells (level 1) or groups (level 2+)
        integer,                                intent(in)  :: nq     ! number of equations for the flow
        real(p2), dimension(nq,ncells),         intent(in)  :: phi    ! dependent variables (du for level 1, del for level 2+)
        real(p2), dimension(:,:,:),             intent(in)  :: V   ! Values (nq x nq block matrix) plus corresponding index
        integer, dimension(:),                  intent(in)  :: C   ! Column index of each value
        integer,dimension(ncells+1),            intent(in)  :: R   ! Start index of each new row
        integer,                                intent(in)  :: nnz
        real(p2), dimension(nq,ncells),         intent(in)  :: res    ! RHS of the NL equation
        integer,                                intent(in)  :: level  ! Multigrid level (not sure if we'll need it yet)
        ! OUTPUT
        integer,                                intent(out) :: ngroup
        integer,                                intent(out) :: nnz_restrict
        integer,  dimension(ncells),            intent(out) :: prolongC
        real(p2), dimension(:,:,:), pointer,    intent(out) :: RAP_V
        integer, dimension(:), pointer,         intent(out)  :: RAP_C
        integer, dimension(:), pointer,         intent(out)  :: RAP_R
        real(p2), dimension(:,:,:), pointer,    intent(out)  :: RAP_Dinv
        real(p2), dimension(:,:), pointer,      intent(out) :: defect_res ! R(A*phi + b)

        ! ! Call the linear solver to relax the newly coarsened system and apply the correction.
        ! allocate(g_correction(5,ngroup))
        ! call linear_sweeps(ngroup, RAP_V,RAP_C,RAP_R,nnz,defect_res,RAP_Dinv, level + 1, g_correction)




        ! real(p2), dimension(5,ncells),          intent(out) :: correction_del ! prolongated correction factor (del) to be 
                                                                              ! passed to finer level

        ! Local Vars
        real(p2), dimension(:,:), allocatable :: g_correction ! restricted correction factor
        real(p2)                            :: strongest_A11, group_size, current_n
        real(p2), dimension(3)              :: A11_nghbrs
        integer,  dimension(3)              :: strongest_3nghbr, sort_index
        integer, dimension(ncells)          :: assign_group
        integer                             :: i, j, k, ck, gi, kk, os, ii, nz_count
        integer                             :: idestat

        real(p2), dimension(nq,ncells)       :: defect ! defect used in RHS for AMG

        ! Restricted Vars
        integer, dimension(ncells)              :: RestrictC ! ProlongC defined above
        integer, dimension(ncells + 1)          :: RestrictR, ProlongR ! Restrict array is longer than needed but this proves fine.

        

        ! Initialize var
        assign_group = 0
        ! correction_del = zero
        ngroup = 0
        RestrictR(1) = 1
        ProlongR = (/ (i, i = 1,(ncells+1)) /)
        assign_cells : do i = 1,ncells
            if ((assign_group(i) /= 0)) cycle
            ! Initialize
            ngroup = ngroup + 1
            group_size = 1
            A11_nghbrs = zero
            strongest_3nghbr = 0
            assign_group(i) = ngroup
            do j = R(i),(R(i+1)-1)
                if (C(j) == i .or. assign_group(C(j)) /= 0) cycle ! skip diagonal and already assigned
                current_n = ABS(V(1,1,j))
                if (current_n > A11_nghbrs(1)) then
                    A11_nghbrs(1) = current_n
                    strongest_3nghbr(1) = C(j)
                    call insertion_sort_real(A11_nghbrs,sort_index)
                    A11_nghbrs(:) = A11_nghbrs(sort_index)
                    strongest_3nghbr(:) = strongest_3nghbr(sort_index)
                end if
            end do
            if ( all(strongest_3nghbr(:) == 0)) then ! all neighbors are in a group
                ! Were just making a group of 1.  May look to add to a neighboring group in the future.  Adds some complexity to the
                ! CSRM foramtion of the restriction matrix.
                ProlongC(i) = ngroup
                RestrictR(ngroup + 1)        = RestrictR(ngroup) + 1 
                RestrictC(RestrictR(ngroup)) = i
                cycle assign_cells
            end if
            ! We only make it here if group_size < 4
            ! So we're just going to keep adding cells till we get to 4
            do j = 1,3
                if ( strongest_3nghbr(j) /= 0 ) then
                    assign_group(strongest_3nghbr(j)) = ngroup
                end if
            end do
            if ( all(strongest_3nghbr(:) /= 0) ) then
                ProlongC(i) = ngroup
                do j = 1,3
                    ProlongC(strongest_3nghbr(j)) = ngroup
                end do
                RestrictR(ngroup + 1) = RestrictR(ngroup) + 4
                RestrictC(RestrictR(ngroup):(RestrictR(ngroup+1)-1)) = insertion_sort_int((/i,strongest_3nghbr/),4)
                cycle ! if all 3 neighbors have been assigned we don't need to look for more
                ! only get here if they haven't been
            end if
            nghbr_group_add : do j = R(strongest_3nghbr(3)) , (R(strongest_3nghbr(3) + 1) - 1) 
                ! move to row of strongest neighbor
                if ( assign_group(C(j)) /= 0 .or. j == C(R(i)) ) cycle ! allready assigned or diagonal
                current_n = ABS(V(1,1,j))
                if (current_n > A11_nghbrs(1)) then
                    A11_nghbrs(1) = current_n
                    strongest_3nghbr(1) = C(j)
                    call insertion_sort_real(A11_nghbrs,sort_index)
                    A11_nghbrs(:) = A11_nghbrs(sort_index)
                    strongest_3nghbr(:) = strongest_3nghbr(sort_index)
                end if
                if ( all(strongest_3nghbr(:) /= 0)) exit nghbr_group_add ! only adding neighbors, not replacing.  Once we fill up
                ! were done.
            end do nghbr_group_add 
            ! add any missing terms
            nz_count = 0
            do j = 1,3
                if ( strongest_3nghbr(j) /= 0 ) then
                    assign_group(strongest_3nghbr(j)) = ngroup
                    ProlongC(strongest_3nghbr(j)) = ngroup
                    nz_count = nz_count + 1
                end if
            end do
            ProlongC(i) = ngroup
            RestrictR(ngroup + 1) = RestrictR(ngroup) + nz_count + 1
            if (nz_count > 0) then
                RestrictC(RestrictR(ngroup):(RestrictR(ngroup+1)-1)) = &
                    insertion_sort_int((/i,strongest_3nghbr((4-nz_count):)/),nz_count + 1)
            else
                RestrictC(RestrictR(ngroup)) = i ! I don't think it's possible to get this far with nz_count = 0 but I'm leaving
                ! this here just in case...
            end if
        end do assign_cells

        ! Create coarse level operetor A^H = RAP
        allocate(RAP_R(ngroup + 1))
        call R_A_P(ncells,ngroup,nq,nnz,RestrictC,RestrictR,ProlongC,ProlongR,V,C,R,RAP_V,RAP_C,RAP_R,nnz_restrict)
        
        ! Create inverse block matrix of RAP diagonal terms
        allocate(RAP_Dinv(5,5,ngroup))
        do gi = 1,ngroup
            do j = RAP_R(gi),(RAP_R(gi+1)-1)
                if (RAP_C(j) == gi) then ! diag
                    call gewp_solve(RAP_V(:,:,j),5, RAP_Dinv(:,:,gi),idestat)
                end if
            end do
        end do
        
        if (idestat/=0) then
            write(*,*) "  Error in inverting the diagonal block... Stop"
            write(*,*) "  Group number = ", i
            write(*,*) "  Level        = ", level
            stop
        endif

        ! Calculate and restrict the defect d = A*phi + b
        call compute_defect(ncells,nq,V,C,R,phi,res,defect)
        allocate(defect_res(5,ngroup))
        defect_res = zero
        do i = 1,ngroup
            do j = RestrictR(i),(RestrictR(i+1)-1)
                defect_res(:,i) = defect_res(:,i) + defect(:,RestrictC(j))
            end do
        end do


        ! ##################################################################
        ! ##################################################################
        ! ##################################################################
        ! ##################################################################

        ! Call the linear solver to relax the newly coarsened system and apply the correction.
        ! allocate(g_correction(5,ngroup))
        ! call linear_sweeps(ngroup, RAP_V,RAP_C,RAP_R,nnz,defect_res,RAP_Dinv, level + 1, g_correction)
        ! do i = 1,ncells
        !     ! Since ProlongC has 1 value per row we can skip the inner j loop.
        !     correction_del(:,i) = g_correction(:,ProlongC(i))
        ! end do

        ! ! Make sure all of the allocated arrays are deallocated
        ! if (allocated(RAP_V))        deallocate(     RAP_V)
        ! if (allocated(RAP_C))        deallocate(     RAP_C)
        ! if (allocated(RAP_R))        deallocate(     RAP_R)
        ! if (allocated(RAP_Dinv))     deallocate(  RAP_Dinv)
        ! if (allocated(defect_res))   deallocate(defect_res)
        ! if (allocated(g_correction)) deallocate(g_correction)

        ! direction = DOWN ! we've finished a multigrid level that means we're going back down the levels
    
    end subroutine algebraic_multigrid_restrict_legacy

    subroutine amg_restric_rs(ncells,nq,phi,V,C,R,nnz,res,level, &
        ngroup,nnz_restrict,prolongC,RAP_V,RAP_C,RAP_R,RAP_Dinv,defect_res)

        use common                  , only : p2, zero
        
        use grid                    , only : cc_data_type
        
        use direct_solve            , only : gewp_solve

        use sparse_common           , only : insertion_sort_int, insertion_sort_real

        use ruge_stuben             , only : rs_agglom


        !use stdlib_sorting, only : sort
        
        implicit none

        ! INPUT
        integer,                                intent(in)  :: ncells ! number of cells (level 1) or groups (level 2+)
        integer,                                intent(in)  :: nq     ! number of equations for the flow
        real(p2), dimension(:,:),               intent(in)  :: phi    ! dependent variables (du for level 1, del for level 2+)
        real(p2), dimension(:,:,:),             intent(in)  :: V   ! Values (nq x nq block matrix) plus corresponding index
        integer, dimension(:),                  intent(in)  :: C   ! Column index of each value
        integer, dimension(:),                  intent(in)  :: R   ! Start index of each new row
        integer,                                intent(in)  :: nnz
        real(p2), dimension(:,:),               intent(in)  :: res    ! RHS of the NL equation
        integer,                                intent(in)  :: level  ! Multigrid level (not sure if we'll need it yet)
        ! OUTPUT
        integer,                                intent(out) :: ngroup
        integer,                                intent(out) :: nnz_restrict
        integer,  dimension(:),                 intent(out) :: prolongC
        real(p2), dimension(:,:,:), pointer,    intent(out) :: RAP_V
        integer, dimension(:), pointer,         intent(out)  :: RAP_C
        integer, dimension(:), pointer,         intent(out)  :: RAP_R
        real(p2), dimension(:,:,:), pointer,    intent(out)  :: RAP_Dinv
        real(p2), dimension(:,:), pointer,      intent(out) :: defect_res ! R(A*phi + b)

        ! Local Vars
        integer                                 :: i, j, gi
        integer                                 :: idestat
        real(p2), dimension(nq,ncells)          :: defect ! defect used in RHS for AMG
        integer, dimension(ncells)              :: RestrictC ! ProlongC defined above
        integer, dimension(ncells + 1)          :: ProlongR ! Restrict array is longer than needed but this proves fine.
        integer, dimension(:), pointer          :: RestrictR


        
        nullify(RestrictR)

        call rs_agglom(ncells, C, R, ngroup, RestrictR, RestrictC, ProlongR, prolongC)

        ! Create coarse level operetor A^H = RAP
        allocate(RAP_R(ngroup + 1))
        call R_A_P(ncells,ngroup,nq,nnz,RestrictC,RestrictR,ProlongC,ProlongR,V,C,R,RAP_V,RAP_C,RAP_R,nnz_restrict)
        
        ! Create inverse block matrix of RAP diagonal terms
        allocate(RAP_Dinv(5,5,ngroup))
        do gi = 1,ngroup
            do j = RAP_R(gi),(RAP_R(gi+1)-1)
                if (RAP_C(j) == gi) then ! diag
                    call gewp_solve(RAP_V(:,:,j),5, RAP_Dinv(:,:,gi),idestat)
                    if (idestat/=0) then
                        write(*,*) " amg_restrict_rs: Error in inverting the diagonal block... Stop"
                        write(*,*) "  Group number = ", i
                        write(*,*) "  Level        = ", level
                        stop
                    endif
                end if
            end do
        end do
        
        

        ! Calculate and restrict the defect d = A*phi + b
        call compute_defect(ncells,nq,V,C,R,phi,res,defect)
        allocate(defect_res(5,ngroup))
        defect_res = zero
        do i = 1,ngroup
            do j = RestrictR(i),(RestrictR(i+1)-1)
                defect_res(:,i) = defect_res(:,i) + defect(:,RestrictC(j))
            end do
        end do

    end subroutine amg_restric_rs


    subroutine algebraic_multigrid_prolong(ncells,prolongC,restricted_correction,correction)
    
        use common , only : p2
        
        implicit none
    
        ! INPUT
        integer,                 intent(in) :: ncells
        integer, dimension(:)  , intent(in) :: prolongC
        real(p2),dimension(:,:), intent(in) :: restricted_correction
        ! OUTPUT
        real(p2),dimension(:,:), intent(inout) :: correction

        integer :: i

        do i = 1,ncells
            ! Since ProlongC has 1 value per row we can skip the inner j loop.
            correction(:,i) = correction(:,i) + restricted_correction(:,ProlongC(i))
        end do

    end subroutine algebraic_multigrid_prolong

    subroutine A_times_P(nq,ncells,nnz,V,C,R,ProlongC,ProlongR,productV,productC,productR)

        use common, only : p2, zero

        use sparse_matrix , only : sparse_sparse_pre_alloc, sparse_real_times_sparse_bool
        ! This subroutine performes the first half (AP) of RAP (where P = R^T).  The preallocation and computation are performed
        ! using the row-wise product.  This takes advantage of the sparse nature of the cells for a very efficient matrix-matrix
        ! product.  Additionally it allows for multiplication without needing to match indices.

        implicit none

        integer, intent(in) :: nq
        integer, intent(in) :: ncells
        integer, intent(in) :: nnz
        real(p2), dimension(nq,nq,nnz), intent(in) :: V
        integer, dimension(nnz), intent(in) :: C
        integer, dimension(ncells+1), intent(in) :: R
        integer, dimension(ncells), intent(in) :: ProlongC
        integer, dimension(ncells + 1), intent(in) :: ProlongR
  

        real(p2), dimension(:,:,:), allocatable, intent(out) :: productV
        integer, dimension(:), allocatable, INTENT(OUT) :: productC
        integer, dimension(ncells + 1), INTENT(OUT) :: productR

        integer :: nnz_prime, i, j, k, os, counter, nnz_prime_new
        logical :: added
        real(p2), dimension(5,5) :: intermediateProduct

        ! Compute the resultant number of nonzero values.
        call sparse_sparse_pre_alloc(ncells,C,R,ProlongC,ProlongR,nnz_prime)

        ! Compute A*Prolong
        allocate(productV(5,5,nnz_prime))
        allocate(productC(    nnz_prime))
        call sparse_real_times_sparse_bool(nq,ncells,V,C,R,ProlongC,ProlongR,nnz_prime,productV,productC,productR)

    end subroutine A_times_P


    subroutine R_A_P(ncells,ngroups,nq,nnz,RestrictC,RestrictR,ProlongC,ProlongR,V,C,R,RAP_V,RAP_C,RAP_R,nnz_prime_final)
        ! This subroutine computes the restriction matrix A^H = RAP for algebraic multi grid and stores it using BCSM.
        use common , only : p2, zero

        use sparse_matrix , only : sparse_sparse_pre_alloc, sparse_bool_times_sparse_real

        implicit none
        integer, intent(in) :: ncells                                   ! # of cells for the coarse mesh
        integer, intent(in) :: ngroups                                  ! # of groups on the restricted level
        integer, intent(in) :: nq                                       ! # of equations in block matrix
        integer, intent(in) :: nnz                                      ! # of nonzero entries in the fine A matrix
        integer, dimension(ncells), intent(in) :: RestrictC             ! Restriction matrix Columns
        integer, dimension(ngroups + 1), intent(in) :: RestrictR        ! Restriction matrix Columns
        integer, dimension(ncells), intent(in) :: ProlongC              ! Restriction matrix Columns
        integer, dimension(ncells + 1), intent(in) :: ProlongR          ! Restriction matrix Columns
        real(p2), dimension(nq,nq,nnz), intent(in) :: V                   ! Block Sparse Compressed Matrix (BSCM) values of A matrix
        integer, dimension(nnz), intent(in) :: C                        ! BCSM Columns of A matrix
        integer, dimension(ncells+1), intent(in) :: R                   ! BCSM Rows of A matrix
        ! integer,dimension(ngroups,ncells) :: Restrict

        real(p2), dimension(:,:,:), pointer, INTENT(OUT) :: RAP_V   ! BCSM Values of restricted A matrix
        integer, dimension(:), pointer, INTENT(OUT) :: RAP_C        ! BCSM Columns of restricted A matrix
        integer, dimension(ngroups + 1), INTENT(OUT) :: RAP_R           ! BCSM Rows of restricted A matrix
        integer, optional, intent(out) :: nnz_prime_final               ! (Optional) # of nonzero values in restricted A matrix

        real(p2), dimension(:,:,:), allocatable :: AP_V                 ! BCSM Values of intermediate A*(R^T) matrix
        integer, dimension(:), allocatable :: AP_C                      ! BCSM Rows of intermediate A*(R^T) matrix
        integer, dimension(ncells + 1) :: AP_R                          ! BCSM Columns of intermediate A*(R^T) matrix

        integer :: i, j, k, nnz_prime

        nnz_prime = 0
        
        call A_times_P(nq,ncells,nnz,V,C,R,ProlongC,ProlongR,AP_V,AP_C,AP_R)

        call sparse_sparse_pre_alloc(ngroups,RestrictC,RestrictR,AP_C,AP_R,nnz_prime)

        allocate(RAP_V(5,5,nnz_prime))
        allocate(RAP_C(    nnz_prime))
        call sparse_bool_times_sparse_real(nq,ngroups,RestrictC,RestrictR,AP_V,AP_C,AP_R,nnz_prime,RAP_V,RAP_C,RAP_R)

        deallocate(AP_V)
        deallocate(AP_C)

        if (present(nnz_prime_final)) nnz_prime_final = nnz_prime

    end subroutine R_A_P 

    subroutine compute_defect(ncells,nq,V,C,R,phi,b,defect)
        ! Computes the defect d = A*phi + b
        use common        , only : p2
        use sparse_matrix , only : sparseblock_times_vectorblock

        implicit none

        integer, intent(in) :: ncells
        integer, intent(in) :: nq
        real(p2), dimension(:,:,:), intent(in) :: V
        integer, dimension(:), intent(in) :: C
        integer, dimension(ncells + 1) :: R
        real(p2), dimension(nq,ncells), intent(in) :: phi
        real(p2), dimension(nq,ncells), intent(in) :: b

        real(p2), dimension(nq,ncells), intent(out) :: defect

        real(p2), dimension(nq,ncells) :: product

        call sparseblock_times_vectorblock(ncells,nq,V,C,R,phi,product)

        defect = product + b

    end subroutine compute_defect

    function convert_amg_c_to_i(amg_char) result(amg_int)

        character(1), intent(in) :: amg_char
        integer                  :: amg_int

        select case(amg_char)
        case('f')
            amg_int = AMG_F
        case('w')
            amg_int = AMG_W
        case('v')
            amg_int = AMG_V
        case default
            write(*,*) "convert_amg_c_to_i: invalid AMG cycle tpye. STOP!"
            stop
        end select

    end function convert_amg_c_to_i


end module algebraic_multigird