module sparse_matrix
    ! This is a module for building and manipulating compressed sparse matrices specific to the linear solver 
    
    implicit none

    public :: build_A_BCSM ! builds A matrix in a Block Compressed Sparse Matrix (BCSM) format

contains
    
    subroutine build_A_BCSM(ncells,cell,jac,V,C,R,nnz)
        ! Takes in a cell structure and a jacobian structure and creates a corresponding A matrix using the Yale meethod:
        ! https://en.wikipedia.org/wiki/Sparse_matrix


        use common      , only : p2
        use grid        , only : cc_data_type
        use jacobian    , only : jacobian_type
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

    subroutine A_times_P(ncells,ngroups,nnz,V,C,R,ProlongC,ProlongR,productV,productC,productR)
        use common, only : p2, zero
        ! This subroutine performes the first half (AP) of RAP (where P = R^T).  The preallocation and computation are performed
        ! using the row-wise product.  This takes advantage of the sparse nature of the cells for a very efficient matrix-matrix
        ! product.  Additionally it allows for multiplication without needing to match indices.

        implicit none

        integer, intent(in) :: ncells
        integer, intent(in) :: ngroups
        integer, intent(in) :: nnz
        real(p2), dimension(5,5,nnz), intent(in) :: V
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
        call sparse_sparse_row_wise_product_AP(ncells,V,C,R,ProlongC,ProlongR,nnz_prime,productV,productC,productR)

    end subroutine A_times_P


    subroutine R_A_P(ncells,ngroups,nnz,RestrictC,RestrictR,ProlongC,ProlongR,V,C,R,RAP_V,RAP_C,RAP_R,nnz_prime_final)
        ! This subroutine computes the restriction matrix A^H = RAP for algebraic multi grid and stores it using BCSM.
        use common , only : p2, zero

        implicit none
        integer, intent(in) :: ncells                                   ! # of cells for the coarse mesh
        integer, intent(in) :: ngroups                                  ! # of groups on the restricted level
        integer, intent(in) :: nnz                                      ! # of nonzero entries in the fine A matrix
        integer, dimension(ncells), intent(in) :: RestrictC             ! Restriction matrix Columns
        integer, dimension(ngroups + 1), intent(in) :: RestrictR        ! Restriction matrix Columns
        integer, dimension(ncells), intent(in) :: ProlongC              ! Restriction matrix Columns
        integer, dimension(ncells + 1), intent(in) :: ProlongR          ! Restriction matrix Columns
        real(p2), dimension(5,5,nnz), intent(in) :: V                   ! Block Sparse Compressed Matrix (BSCM) values of A matrix
        integer, dimension(nnz), intent(in) :: C                        ! BCSM Columns of A matrix
        integer, dimension(ncells+1), intent(in) :: R                   ! BCSM Rows of A matrix
        ! integer,dimension(ngroups,ncells) :: Restrict

        real(p2), dimension(:,:,:), allocatable, INTENT(OUT) :: RAP_V   ! BCSM Values of restricted A matrix
        integer, dimension(:), allocatable, INTENT(OUT) :: RAP_C        ! BCSM Columns of restricted A matrix
        integer, dimension(ngroups + 1), INTENT(OUT) :: RAP_R           ! BCSM Rows of restricted A matrix
        integer, optional, intent(out) :: nnz_prime_final               ! (Optional) # of nonzero values in restricted A matrix

        real(p2), dimension(:,:,:), allocatable :: AP_V                 ! BCSM Values of intermediate A*(R^T) matrix
        integer, dimension(:), allocatable :: AP_C                      ! BCSM Rows of intermediate A*(R^T) matrix
        integer, dimension(ncells + 1) :: AP_R                          ! BCSM Columns of intermediate A*(R^T) matrix

        integer :: i, j, k, nnz_prime

        nnz_prime = 0
        
        call A_times_P(ncells,ngroups,nnz,V,C,R,ProlongC,ProlongR,AP_V,AP_C,AP_R)

        call sparse_sparse_pre_alloc(ngroups,RestrictC,RestrictR,AP_C,AP_R,nnz_prime)

        allocate(RAP_V(5,5,nnz_prime))
        allocate(RAP_C(    nnz_prime))
        call sparse_sparse_row_wise_product_RAP(ngroups,RestrictC,RestrictR,AP_V,AP_C,AP_R,nnz_prime,RAP_V,RAP_C,RAP_R)

        deallocate(AP_V)
        deallocate(AP_C)

        if (present(nnz_prime_final)) nnz_prime_final = nnz_prime

    end subroutine R_A_P 

    subroutine sparse_sparse_pre_alloc(A_m,AC,AR,BC,BR,nnz_prime)
        ! This subroutine uses a row-wise product to determine the number of nonzero terms that will result from the multiplication
        ! of A*B.  Because we are only concered in the number of values we do not actually need to evaluate the final value.  
        ! Instead we only need to check if there is a nonzero value at a given index.

        use sparse_common , only : queue_sort
        implicit none
        
        integer, intent(in) :: A_m  ! Height of A
        integer, dimension(:), intent(in) :: AC         ! Column indices of A
        integer, dimension(A_m+1), intent(in) :: AR     ! Row indices of A
    
        integer, dimension(:), intent(in) :: BC         ! Column indices of B
        integer, dimension(:), intent(in) :: BR         ! Row indices of B

        integer, intent(out) :: nnz_prime   ! number of nonzero values in the output

        integer :: i,j,b_i,b_j, column
        integer, dimension(:), pointer          :: queue1_C
        integer, dimension(:), allocatable      :: queue2_C
        allocate(queue1_C(10))
        allocate(queue2_C(4))
        nnz_prime = 0
        do i = 1,A_m
            do j = AR(i),(AR(i+1)-1)
                b_i = AC(j);
                column = 1
                if (j == AR(i)) then
                    if (BR(b_i+1)-BR(b_i) > size(queue1_C)) then
                        if (associated(queue1_C)) deallocate(queue1_C)
                        allocate(queue1_C(BR(b_i+1)-BR(b_i) ))
                    end if
                    queue1_C = 0
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue1_C(column) = BC(b_j)
                        column = column + 1
                    end do
                else
                    if (BR(b_i+1)-BR(b_i) > size(queue2_C)) then
                        if (allocated(queue2_C)) deallocate(queue2_C)
                        allocate(queue2_C(BR(b_i+1)-BR(b_i) ))
                    end if
                    queue2_C = 0
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue2_C(column) = BC(b_j)
                        column = column + 1
                    end do
                    ! Merge the two queues
                    call queue_sort(queue1_C,size(queue1_C),queue2_C,size(queue2_C))
                end if
            end do
            if ((AR(i+1)-AR(i)) == 1) then
                nnz_prime = nnz_prime + column - 1
            else
                nnz_prime = nnz_prime + size(queue1_C)
            end if
        end do

        if (associated(queue1_C)) deallocate(queue1_C)
        if (allocated(queue2_C))  deallocate(queue2_C)
    end subroutine sparse_sparse_pre_alloc

    subroutine sparse_sparse_row_wise_product_RAP(A_m,AC,AR,BV,BC,BR,nnz_prime, prodV,prodC,prodR)
        ! This function performs an efficient computation of A*B=C where A is m*n block matrix and B is n*k logical array, and A,B, 
        ! and C are all stored using a Compressed Sparse Row Matrix format.  The method is based off the paper 
        ! https://doi.org/10.1109/MICRO50266.2020.00068
        ! 
        !
        ! Each row of C is defined as follows:
        ! C(1,:) = A(1,1)*B(1,:) + A(1,2)*B(2,:) + ... + A(1,n)*B(n,:)
        ! This has the benefit of being able to use two CSR matrices without transposition.  Additionally, for sparse matrices the
        ! number of terms in each row is relatively small.  As a result, only nonzero products will be calculated.
        ! It does require merge-soring the resulting rows as they are computed, however because the two queues are already 
        ! pre-sorted this can be done rather efficiently.
        ! 
        ! Also note: Because A is a logical array, and therefore only true indices are stored, the AV vector is unnecessary.
        
        use common , only : p2, zero

        implicit none

        integer, intent(in) :: A_m  ! Height of A
        ! logical, dimension(:), intent(in) :: AV         ! values of A
        integer, dimension(:), intent(in) :: AC         ! Column indices of A
        integer, dimension(A_m+1), intent(in) :: AR     ! Row indices of A
        
        real(p2), dimension(:,:,:), intent(in) :: BV    ! Values of B
        integer, dimension(:), intent(in) :: BC         ! Column indices of B
        integer, dimension(:), intent(in) :: BR         ! Row indices of B

        integer, intent(in) :: nnz_prime    ! number of nonzero values in the output

        real(p2),dimension(5,5,nnz_prime), intent(out) :: prodV ! Values of C
        integer, dimension(nnz_prime), INTENT(OUT) :: prodC     ! Column indices of C
        integer,dimension(A_m+1), INTENT(OUT) :: prodR          ! Row indices of C

        integer :: i,j,b_i,b_j, column, length
        logical :: queue1_link = .false.
        real(p2), dimension(:,:,:), pointer     :: queue1_V
        real(p2), dimension(:,:,:), allocatable :: queue2_V
        integer, dimension(:), pointer          :: queue1_C
        integer, dimension(:), allocatable      :: queue2_C

        integer, dimension(:), allocatable, target      :: tempC
        real(p2), dimension(:,:,:), allocatable, target :: tempV
        
        prodR(1) = 1
        allocate(queue1_C(1))
        allocate(queue2_C(10))
        allocate(queue1_V(5,5,1))
        allocate(queue2_V(5,5,10))
        queue1_C = 0
        queue2_C = 0
        queue1_V = zero
        queue2_V = zero
        do i = 1,A_m
            do j = AR(i),(AR(i+1)-1)
                b_i = AC(j);
                column = 1
                if (j == AR(i)) then
                    if (BR(b_i+1)-BR(b_i) > size(queue1_C)) then
                        if (associated(queue1_C)) deallocate(queue1_C)
                        if (associated(queue1_V)) deallocate(queue1_V)
                        allocate(queue1_C(BR(b_i+1)-BR(b_i) ))
                        allocate(queue1_V(5,5,BR(b_i+1)-BR(b_i) ))
                    end if
                    queue1_C = 0
                    queue1_V = zero
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue1_C(column) = BC(b_j)
                        queue1_V(:,:,column) = BV(:,:,b_j)
                        column = column + 1
                    end do
                else
                    if (BR(b_i+1)-BR(b_i) > size(queue2_C)) then
                        if (allocated(queue2_C)) deallocate(queue2_C)
                        if (allocated(queue2_V)) deallocate(queue2_V)
                        allocate(queue2_C(BR(b_i+1)-BR(b_i) ))
                        allocate(queue2_V(5,5,BR(b_i+1)-BR(b_i) ))
                    end if
                    queue2_C = 0
                    queue2_V = zero
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue2_C(column) = BC(b_j)
                        queue2_V(:,:,column) = BV(:,:,b_j)
                        column = column + 1
                    end do
                    call queue_sort_5x5(queue1_V,queue1_C,size(queue1_C),queue2_V,queue2_C,size(queue2_C))
                end if
            end do
            if ((AR(i+1)-AR(i)) == 1) then
                prodR(i+1) = prodR(i) + column - 1
                prodV(:,:,prodR(i):(prodR(i+1)-1)) = queue1_V(:,:,1:(column - 1))
                prodC(prodR(i):(prodR(i+1)-1)) = queue1_C(1:(column - 1))
            else
                prodR(i+1) = prodR(i) + size(queue1_C)
                prodV(:,:,prodR(i):(prodR(i+1)-1)) = queue1_V
                prodC(prodR(i):(prodR(i+1)-1)) = queue1_C
            end if
        end do

        if (associated(queue1_C)) deallocate(queue1_C)
        if (associated(queue1_V)) deallocate(queue1_V)
        if (allocated(queue2_C))  deallocate(queue2_C)
        if (allocated(queue2_V))  deallocate(queue2_V)
    end subroutine sparse_sparse_row_wise_product_RAP

    subroutine sparse_sparse_row_wise_product_AP(A_m,AV,AC,AR,BC,BR,nnz_prime, prodV,prodC,prodR)
        ! This function performs an efficient computation of A*B=C where A is m*n block matrix and B is n*k logical array, and A,B, 
        ! and C are all stored using a Compressed Sparse Row Matrix format.  The method is based off the paper 
        ! https://doi.org/10.1109/MICRO50266.2020.00068
        ! 
        !
        ! Each row of C is defined as follows:
        ! C(1,:) = A(1,1)*B(1,:) + A(1,2)*B(2,:) + ... + A(1,n)*B(n,:)
        ! This has the benefit of being able to use two CSR matrices without transposition.  Additionally, for sparse matrices the
        ! number of terms in each row is relatively small.  As a result, only nonzero products will be calculated.
        ! It does require merge-soring the resulting rows as they are computed, however because the two queues are already 
        ! pre-sorted this can be done rather efficiently.
        ! 
        ! Also note: Because B is a logical array, and therefore only true indices are stored, the BV vector is unnecessary.
        
        use common , only : p2, zero

        implicit none

        integer, intent(in) :: A_m  ! Height of A
        real(p2), dimension(:,:,:), intent(in) :: AV    ! values of A
        integer, dimension(:), intent(in) :: AC         ! Column indices of A
        integer, dimension(A_m+1), intent(in) :: AR     ! Row indices of A
        
        ! logical, dimension(:), intent(in) :: BV    ! Values of B
        integer, dimension(:), intent(in) :: BC         ! Column indices of B
        integer, dimension(:), intent(in) :: BR         ! Row indices of B

        integer, intent(in) :: nnz_prime    ! number of nonzero values in the output

        real(p2),dimension(5,5,nnz_prime), intent(out) :: prodV ! Values of C
        integer, dimension(nnz_prime), INTENT(OUT) :: prodC     ! Column indices of C
        integer,dimension(A_m+1), INTENT(OUT) :: prodR          ! Row indices of C

        integer :: i,j,b_i,b_j, column
        real(p2), dimension(:,:,:), pointer     :: queue1_V
        real(p2), dimension(:,:,:), allocatable :: queue2_V
        integer, dimension(:), pointer          :: queue1_C
        integer, dimension(:), allocatable      :: queue2_C

        prodR(1) = 1
        allocate(queue1_C(10))
        allocate(queue2_C(10))
        allocate(queue1_V(5,5,10))
        allocate(queue2_V(5,5,10))
        queue1_C = 0
        queue2_C = 0
        queue1_V = zero
        queue2_V = zero
        do i = 1,A_m
            do j = AR(i),(AR(i+1)-1)
                b_i = AC(j);
                column = 1
                if (j == AR(i)) then
                    if (BR(b_i+1)-BR(b_i) > size(queue1_C)) then
                        if (associated(queue1_C)) deallocate(queue1_C)
                        if (associated(queue1_V)) deallocate(queue1_V)
                        allocate(queue1_C(BR(b_i+1)-BR(b_i) ))
                        allocate(queue1_V(5,5,BR(b_i+1)-BR(b_i) ))
                    end if
                    queue1_C = 0
                    queue1_V = zero
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue1_C(column) = BC(b_j)
                        queue1_V(:,:,column) = AV(:,:,j)
                        column = column + 1
                    end do
                else
                    if (BR(b_i+1)-BR(b_i) > size(queue2_C)) then
                        if (allocated(queue2_C)) deallocate(queue2_C)
                        if (allocated(queue2_V)) deallocate(queue2_V)
                        allocate(queue2_C(BR(b_i+1)-BR(b_i) ))
                        allocate(queue2_V(5,5,BR(b_i+1)-BR(b_i) ))
                    end if
                    queue2_C = 0
                    queue2_V = zero
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue2_C(column) = BC(b_j)
                        queue2_V(:,:,column) = AV(:,:,j)
                        column = column + 1
                    end do
                    call queue_sort_5x5(queue1_V,queue1_C,size(queue1_C),queue2_V,queue2_C,size(queue2_C))
                end if
            end do
            if ((AR(i+1)-AR(i)) == 1) then
                prodR(i+1) = prodR(i) + column - 1
                prodV(:,:,prodR(i):(prodR(i+1)-1)) = queue1_V(:,:,1:(column - 1))
                prodC(prodR(i):(prodR(i+1)-1)) = queue1_C(1:(column - 1))
            else
                prodR(i+1) = prodR(i) + size(queue1_C)
                prodV(:,:,prodR(i):(prodR(i+1)-1)) = queue1_V
                prodC(prodR(i):(prodR(i+1)-1)) = queue1_C
            end if
        end do
        if (associated(queue1_C)) deallocate(queue1_C)
        if (associated(queue1_V)) deallocate(queue1_V)
        if (allocated(queue2_C))  deallocate(queue2_C)
        if (allocated(queue2_V))  deallocate(queue2_V)
    end subroutine sparse_sparse_row_wise_product_AP

    subroutine queue_sort_5x5(queue1_V,queue1_C,length1,queue2_V,queue2_C,length2)
        ! This subroutine sorts two pairs of lists, containing a row of values and corresponding column indices, and sorts them 
        ! based on their column indices into a single merged list pair. It takes advantage of the fact both lists are pre sorted. 
        ! Starting fith the first index of each list pair, the column indices are compared.  If they are equaled the values are 
        ! added together and placed in the next availible spot in the merged list with the corresponding column index.  If they are
        ! not equal the value and column index of the smaller column is added to the merged list.
        ! This subroutine only has to go through each list once making it very efficient.
        ! It also has some logic to deal with trailing zeros in either list.
        use common , only : p2, zero
        
        implicit none

        real(p2), dimension(:,:,:), pointer, intent(inout) :: queue1_V
        integer, dimension(:), pointer, intent(inout) :: queue1_C
        integer, intent(in) :: length1
        integer, intent(in) :: length2
        real(p2), dimension(5,5,length2), intent(inout) :: queue2_V
        integer, dimension(length2), intent(inout) :: queue2_C
        

        real(p2), dimension(:,:,:), pointer  :: helperV
        integer, dimension(:), pointer :: helperC
        integer :: i, j, k

        allocate(helperV(5,5,length1+length2))
        allocate(helperC(length1+length2))
        helperV = zero
        helperC = 0
        i = 1
        j = 1
        k = 1
        do
            if (queue1_C(i) == 0 .and. queue2_C(j) == 0) exit
            if (queue1_C(i) == 0) then
                helperV(:,:,k) = queue2_V(:,:,j)
                helperC(k)     = queue2_C(j)
                j = j + 1
                if (j > length2) then
                    j = length2
                    queue2_C(j) = 0
                end if
            elseif (queue2_C(j) == 0) then
                helperV(:,:,k) = queue1_V(:,:,i)
                helperC(k)     = queue1_C(i)
                i = i + 1
                if (i > length1) then
                    i = length1
                    queue1_C(i) = 0
                end if
            elseif (queue1_C(i) > queue2_C(j)) then
                helperV(:,:,k) = queue2_V(:,:,j)
                helperC(k)     = queue2_C(j)
                j = j + 1
                if (j > length2) then
                    j = length2
                    queue2_C(j) = 0
                end if
            elseif (queue1_C(i) == queue2_C(j)) then
                helperV(:,:,k) = queue1_V(:,:,i) + queue2_V(:,:,j)
                helperC(k)     = queue2_C(j)
                i = i + 1
                if (i > length1) then
                    i = length1
                    queue1_C(i) = 0
                end if
                j = j + 1
                if (j > length2) then
                    j = length2
                    queue2_C(j) = 0
                end if
            else
                helperV(:,:,k) = queue1_V(:,:,i)
                helperC(k)     = queue1_C(i)
                i = i + 1
                if (i > length1) then
                    i = length1
                    queue1_C(i) = 0
                end if
            end if
            k = k + 1
        end do
        if (associated(queue1_C)) deallocate(queue1_C)
        if (associated(queue1_V)) deallocate(queue1_V)
        queue1_C => helperC(1:k-1)
        queue1_V => helperV(:,:,1:k-1)
        nullify(helperC)
        nullify(helperV)
        queue2_C = 0
        queue2_V = zero

    end subroutine queue_sort_5x5

    

    subroutine sparseblock_times_vectorblock(ncells,V,C,R,x,b)
        ! This subroutine computes A*x = b where A is a sparse block matrix, represented by V, C, and R, and x and b are dense 
        ! block vectors (a ncells x 1 vector or nEQ x 1 blocks)
        use common , only : p2, zero

        implicit none
        
        integer, intent(in)                    :: ncells
        real(p2), dimension(:,:,:), intent(in) :: V
        integer, dimension(:), intent(in)      :: C
        integer, dimension(ncells+1), intent(in)      :: R
        real(p2), dimension(5,ncells), intent(in)     :: x

        real(p2), dimension(5,ncells), intent(out)    :: b

        real(p2), dimension(5,5)  :: currentV
        real(p2), dimension(5)  :: currentX

        integer :: i,j,ii

        b = zero
        do i = 1,ncells
            do j = R(i),(R(i+1)-1)
                b(:,i) = b(:,i) + matmul(V(:,:,j),x(:,C(j)));
            end do
        end do
        
    end subroutine sparseblock_times_vectorblock

    subroutine vector_times_sparse(ncells,V,C,R,x,b)
        ! ! This subroutine computes x*A = b where A is a sparse block matrix, represented by V, C, and R, and x and b are dense 
        ! block vectors (a 1 x ncells vector or 1 x nEQ blocks)
        use common , only : p2

        implicit none
        
        integer, intent(in)                    :: ncells
        real(p2), dimension(:), intent(in) :: V
        integer, dimension(:), intent(in)      :: C
        integer, dimension(ncells+1), intent(in)      :: R
        real(p2), dimension(ncells), intent(in)     :: x

        real(p2), dimension(ncells), intent(out)    :: b

        integer :: i, j

        do i = 1,ncells
            do j = R(i),(R(i+1)-1)
                b(C(j)) = b(C(j)) +  x(i) * V(j)
            end do
        end do

    end subroutine vector_times_sparse
    

end module sparse_matrix