module sparse_common
    ! a set of subroutines that are common to both block and scalar sparse matrix operations

    implicit none

    contains

    subroutine queue_sort(queue1_C,length1,queue2_C,length2)
        ! This subroutine sorts two lists in the same method as queue_sort_5x5.  The only difference is it only takes in lists of 
        ! column indices (no value lists).
        use common , only : p2, zero
        
        implicit none

        integer, dimension(:), pointer, intent(inout) :: queue1_C
        integer, intent(in) :: length1
        integer, intent(in) :: length2
        integer, dimension(length2), intent(inout) :: queue2_C
        
        integer, dimension(:), pointer :: helperC
        integer :: i, j, k

        allocate(helperC(length1+length2))
        helperC = 0
        i = 1
        j = 1
        k = 1
        do
            if (queue1_C(i) == 0 .and. queue2_C(j) == 0) exit
            if (queue1_C(i) == 0) then
                helperC(k)     = queue2_C(j)
                j = j + 1
                if (j > length2) then
                    j = length2
                    queue2_C(j) = 0
                end if
            elseif (queue2_C(j) == 0) then
                helperC(k)     = queue1_C(i)
                i = i + 1
                if (i > length1) then
                    i = length1
                    queue1_C(i) = 0
                end if
            elseif (queue1_C(i) > queue2_C(j)) then
                helperC(k)     = queue2_C(j)
                j = j + 1
                if (j > length2) then
                    j = length2
                    queue2_C(j) = 0
                end if
            elseif (queue1_C(i) == queue2_C(j)) then
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
        queue1_C => helperC(1:k-1)
        nullify(helperC)
        queue2_C = 0
        
    end subroutine queue_sort

    subroutine vector_times_sparse(ncells,V,C,R,x,b)
        ! ! This subroutine computes x*A = b where A is a sparse block matrix, represented by V, C, and R, and x and b are dense 
        ! block vectors (a 1 x ncells vector or 1 x nEQ blocks)
        use common , only : p2

        implicit none
        
        integer,                intent(in) :: ncells
        real(p2), dimension(:), intent(in) :: V
        integer,  dimension(:), intent(in) :: C
        integer,  dimension(:), intent(in) :: R
        real(p2), dimension(:), intent(in) :: x

        real(p2), dimension(:), intent(out):: b

        integer :: i, j

        do i = 1,ncells
            do j = R(i),(R(i+1)-1)
                b(C(j)) = b(C(j)) +  x(i) * V(j)
            end do
        end do

    end subroutine vector_times_sparse

    subroutine sparse_sparse_pre_alloc(A_m,AC,AR,BC,BR,nnz_prime)
        ! This subroutine uses a row-wise product to determine the number of nonzero terms that will result from the multiplication
        ! of A*B.  Because we are only concered in the number of values we do not actually need to evaluate the final value.  
        ! Instead we only need to check if there is a nonzero value at a given index.

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

    subroutine insertion_sort_index(x,index)
        ! This routine uses insertion sort sort the input vector x and returns the output vector index
        ! which is the index needed to sort x.
        ! Insertion sort is O(n^2) which is generally inefficient.  However its simplicity allows it 
        ! to be faster than most O(n log n) methods for small arrays.  Since the incoming array will 
        ! just be the cell neighbors for the coarse mesh n <= 6 unless I get around to adding support
        ! for general polyhedral cells (please no...)
        !
        ! ref: https://en.wikipedia.org/wiki/Insertion_sort

        integer, dimension(:), intent(in)  :: x
        integer, dimension(:), intent(out) :: index
        
        integer, dimension(:), allocatable :: x_internal
        integer i, j, length, x_j, index_j
        length = size(x)
        allocate(x_internal(length))
        x_internal = x

        index = (/(i,i=1,length)/)
        
        do i = 1,length
            j = i
            inner : do 
                if ( (j <= 1) ) then
                    exit inner
                else if ( x_internal(j-1) < x_internal(j)) then ! it doesn't like evaluating this line with .or. 
                    !if j = 1 (array bounds)
                    ! I know this is stupid, but I don't want to debug it right now...
                    exit inner
                end if
                x_j = x_internal(j)
                x_internal(j) = x_internal(j-1)
                x_internal(j-1) = x_j
                index_j = index(j)
                index(j) = index(j-1)
                index(j-1) = index_j
                j = j-1
            end do inner
        end do
        deallocate(x_internal)
    end subroutine insertion_sort_index

    function insertion_sort_int(x,length)
        ! Performs insertion sort on the list x of integers
        ! While insertion sort is O(n^2) it is generally very efficient for very short lists. The maximum number of elements this
        ! subroutine will experience is currently 4 so a very simple sorting method should be more efficient than a better scaling
        ! but more complex method.
        
        implicit none
        integer, intent(in)                    :: length
        integer, dimension(length), intent(in) :: x
        integer, dimension(length)             :: insertion_sort_int

        integer :: i,j,x_j
        insertion_sort_int = x

        do i = 1,length
            j = i
            inner : do
                if ( (j <= 1) ) then
                    exit inner
                else if (insertion_sort_int(j-1) < insertion_sort_int(j)) then
                    exit inner
                end if
                x_j = insertion_sort_int(j)
                insertion_sort_int(j) = insertion_sort_int(j-1)
                insertion_sort_int(j-1) = x_j
                j = j-1
            end do inner
        end do
    end function

    subroutine insertion_sort_real(x,index)
        use common , only : p2
        ! This routine uses insertion sort sort the input vector x and returns the output vector index
        ! which is the index needed to sort x.
        ! Insertion sort is O(n^2) which is generally inefficient.  However its simplicity allows it 
        ! to be faster than most O(n log n) methods for small arrays.  Since the incoming array will 
        ! just be the cell neighbors for the coarse mesh n <= 6 unless I get around to adding support
        ! for general polyhedral cells (please no...)
        !
        ! ref: https://en.wikipedia.org/wiki/Insertion_sort

        real(p2), dimension(:), intent(in)  :: x
        integer, dimension(:), intent(out) :: index
        
        real(p2), dimension(:), allocatable :: x_internal
        integer i, j, length, index_j
        real(p2) x_j
        length = size(x)
        allocate(x_internal(length))
        x_internal = x

        index = (/(i,i=1,length)/)
        
        do i = 1,length
            j = i
            inner : do 
                if ( (j <= 1) ) then
                    exit inner
                else if ( x_internal(j-1) < x_internal(j)) then ! it doesn't like evaluating this line with .or. if j = 1 (array bounds)
                                              ! I know this is stupid, but I don't want to debug it right now...
                    exit inner
                end if
                x_j = x_internal(j)
                x_internal(j) = x_internal(j-1)
                x_internal(j-1) = x_j
                index_j = index(j)
                index(j) = index(j-1)
                index(j-1) = index_j
                j = j-1
            end do inner
        end do
        deallocate(x_internal)
    end subroutine insertion_sort_real


    subroutine destroy_A(V,C)
        ! This is supposed to be a catch all subroutine for deallocating V and C but for some reason fortran got grumpy when I 
        ! tried to use it so I haven't been...
        use common , only : p2
        ! Deallocates any allocated value and column arrays
        real(p2), dimension(:,:,:),optional, allocatable, intent(out)  :: V   ! Values (5x5 block matrix) plus corresponding index
        integer, dimension(:),optional, allocatable, intent(out)       :: C   ! Column index of each value

        if (allocated(V)) deallocate(V)
        if (allocated(C)) deallocate(C)

    end subroutine
end module sparse_common