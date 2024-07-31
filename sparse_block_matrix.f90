module sparse_matrix
    ! This is a module for building and manipulating compressed sparse matrices specific to the linear solver 
    
    implicit none

    public

    ! This will allow us to abstract the scalar and block operations
    interface sparse_bool_times_sparse_real
        module procedure sparse_bool_times_sparse_real_block
    end interface sparse_bool_times_sparse_real

    interface sparse_real_times_sparse_bool
        module procedure sparse_real_times_sparse_bool_block
    end interface sparse_real_times_sparse_bool
contains
    
    

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

    subroutine sparse_bool_times_sparse_real_block(nq,A_m,AC,AR,BV,BC,BR,nnz_prime, prodV,prodC,prodR)
        ! This function performs an efficient computation of A*B=C where A is m*n block matrix and B is n*k logical array, and A,B, 
        ! and C are all stored using a Compressed Sparse Row Matrix format.  The method is based off the paper 
        ! https://doi.org/10.1109/MICRO50266.2020.00068
        ! 
        ! NOTE: The vector A does not actually contain any logical variables.  The vector of logical values *would* be AV; however,
        ! because AV only stores the nonzero array values every value would be true.  Therefore we only need the indices of the 
        ! nonzero terms (AC and AR) to form A
        !
        !
        ! Each row of C is defined as follows:
        ! C(1,:) = A(1,1)*B(1,:) + A(1,2)*B(2,:) + ... + A(1,n)*B(n,:)
        ! This has the benefit of being able to use two CSR matrices without transposition.  Additionally, for sparse matrices the
        ! number of terms in each row is relatively small.  As a result, only nonzero products will be calculated.
        ! It does require merge-soring the resulting rows as they are computed, however because the two queues are already 
        ! pre-sorted this can be done rather efficiently.
        ! 
        ! This array is used for computing by the algebraic multigrid module to multiply R*(AP) where AP is the (n x m) Jacobian  
        ! and has already been multiplied by P the (n x m) prolongation matrix equal to the transpose of the Restriction matrix R.
        
        use common , only : p2, zero

        implicit none

        ! INPUT
        integer, intent(in) :: nq
        integer, intent(in) :: A_m  ! Height of A
        ! logical, dimension(:), intent(in) :: AV         ! values of A
        integer, dimension(:), intent(in) :: AC         ! Column indices of A
        integer, dimension(A_m+1), intent(in) :: AR     ! Row indices of A
        
        real(p2), dimension(:,:,:), intent(in) :: BV    ! Values of B
        integer, dimension(:), intent(in) :: BC         ! Column indices of B
        integer, dimension(:), intent(in) :: BR         ! Row indices of B

        integer, intent(in) :: nnz_prime    ! number of nonzero values in the output

        real(p2),dimension(nq,nq,nnz_prime), intent(out) :: prodV ! Values of C
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
        allocate(queue1_V(nq,nq,1))
        allocate(queue2_V(nq,nq,10))
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
                        allocate(queue1_V(nq,nq,BR(b_i+1)-BR(b_i) ))
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
                        allocate(queue2_V(nq,nq,BR(b_i+1)-BR(b_i) ))
                    end if
                    queue2_C = 0
                    queue2_V = zero
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue2_C(column) = BC(b_j)
                        queue2_V(:,:,column) = BV(:,:,b_j)
                        column = column + 1
                    end do
                    call queue_sort_block(nq,queue1_V,queue1_C,size(queue1_C),queue2_V,queue2_C,size(queue2_C))
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
    end subroutine sparse_bool_times_sparse_real_block

    subroutine sparse_real_times_sparse_bool_block(nq,A_m,AV,AC,AR,BC,BR,nnz_prime, prodV,prodC,prodR)
        ! This function performs an efficient computation of A*B=C where A is m*n block matrix and B is n*k logical array, and A,B, 
        ! and C are all stored using a Compressed Sparse Row Matrix format.  The method is based off the paper 
        ! https://doi.org/10.1109/MICRO50266.2020.00068
        ! 
        ! NOTE: The vector B does not actually contain any logical variables.  The vector of logical values *would* be BV; however,
        ! because BV only stores the nonzero array values every value would be true.  Therefore we only need the indices of the 
        ! nonzero terms (BC and BR) to form B
        !
        ! Each row of C is defined as follows:
        ! C(1,:) = A(1,1)*B(1,:) + A(1,2)*B(2,:) + ... + A(1,n)*B(n,:)
        ! This has the benefit of being able to use two CSR matrices without transposition.  Additionally, for sparse matrices the
        ! number of terms in each row is relatively small.  As a result, only nonzero products will be calculated.
        ! It does require merge-soring the resulting rows as they are computed, however because the two queues are already 
        ! pre-sorted this can be done rather efficiently.
        ! 
        ! This array is used for computing by the algebraic multigrid module to multiply A*P where A is the (n x n) Jacobian and P 
        ! is the (n x m) prolongation matrix equal to the transpose of the Restriction matrix R.
        
        use common , only : p2, zero

        implicit none

        integer, intent(in) :: nq
        integer, intent(in) :: A_m  ! Height of A
        real(p2), dimension(:,:,:), intent(in) :: AV    ! values of A
        integer, dimension(:), intent(in) :: AC         ! Column indices of A
        integer, dimension(A_m+1), intent(in) :: AR     ! Row indices of A
        
        ! logical, dimension(:), intent(in) :: BV    ! Values of B
        integer, dimension(:), intent(in) :: BC         ! Column indices of B
        integer, dimension(:), intent(in) :: BR         ! Row indices of B

        integer, intent(in) :: nnz_prime    ! number of nonzero values in the output

        real(p2),dimension(nq,nq,nnz_prime), intent(out) :: prodV ! Values of C
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
        allocate(queue1_V(nq,nq,10))
        allocate(queue2_V(nq,nq,10))
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
                        allocate(queue1_V(nq,nq,BR(b_i+1)-BR(b_i) ))
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
                        allocate(queue2_V(nq,nq,BR(b_i+1)-BR(b_i) ))
                    end if
                    queue2_C = 0
                    queue2_V = zero
                    do b_j = BR(b_i),(BR(b_i+1)-1)
                        queue2_C(column) = BC(b_j)
                        queue2_V(:,:,column) = AV(:,:,j)
                        column = column + 1
                    end do
                    call queue_sort_block(nq,queue1_V,queue1_C,size(queue1_C),queue2_V,queue2_C,size(queue2_C))
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
    end subroutine sparse_real_times_sparse_bool_block

    subroutine sparseblock_times_vectorblock(ncells,nq,V,C,R,x,b)
        ! This subroutine computes A*x = b where A is a sparse block matrix, represented by V, C, and R, and x and b are dense 
        ! block vectors (a ncells x 1 vector or nEQ x 1 blocks)
        use common , only : p2, zero

        implicit none
        
        integer, intent(in)                    :: ncells
        integer, intent(in)                    :: nq
        real(p2), dimension(:,:,:), intent(in) :: V
        integer, dimension(:), intent(in)      :: C
        integer, dimension(ncells+1), intent(in)      :: R
        real(p2), dimension(nq,ncells), intent(in)     :: x

        real(p2), dimension(nq,ncells), intent(out)    :: b

        real(p2), dimension(nq,nq)  :: currentV
        real(p2), dimension(nq)  :: currentX

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
    
    subroutine queue_sort_block(nq,queue1_V,queue1_C,length1,queue2_V,queue2_C,length2)
        ! This subroutine sorts two pairs of lists, containing a row of values and corresponding column indices, and sorts them 
        ! based on their column indices into a single merged list pair. It takes advantage of the fact both lists are pre sorted. 
        ! Starting fith the first index of each list pair, the column indices are compared.  If they are equaled the values are 
        ! added together and placed in the next availible spot in the merged list with the corresponding column index.  If they are
        ! not equal the value and column index of the smaller column is added to the merged list.
        ! This subroutine only has to go through each list once making it very efficient.
        ! It also has some logic to deal with trailing zeros in either list.
        use common , only : p2, zero
        
        implicit none

        integer, intent(in) :: nq
        real(p2), dimension(:,:,:), pointer, intent(inout) :: queue1_V
        integer, dimension(:), pointer, intent(inout) :: queue1_C
        integer, intent(in) :: length1
        integer, intent(in) :: length2
        real(p2), dimension(nq,nq,length2), intent(inout) :: queue2_V
        integer, dimension(length2), intent(inout) :: queue2_C
        

        real(p2), dimension(:,:,:), pointer  :: helperV
        integer, dimension(:), pointer :: helperC
        integer :: i, j, k

        allocate(helperV(nq,nq,length1+length2))
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

    end subroutine queue_sort_block

end module sparse_matrix