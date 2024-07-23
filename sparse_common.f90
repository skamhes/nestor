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