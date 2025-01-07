module sorting

    use common , only : p2

    implicit none

    private

    public heap_sort_index

    contains

    subroutine heap_sort_index(values, index, count, sorted_index)

        ! based off the pseudo-code from https://en.wikipedia.org/wiki/Heapsort#Standard_implementation
        ! note: because the code from wiki uses 0 indexing, some of the integers differ by 1 in this implementation.
        implicit none

        real(p2), dimension(:), intent(in)  :: values
        integer,  dimension(:), intent(in)  :: index
        integer,                intent(in)  :: count
        integer,  dimension(:), intent(out) :: sorted_index

        real(p2), dimension(count) :: sorted_values

        integer :: istart, iend, iroot, ichild

        sorted_index  = index
        sorted_values = values

        istart = (count/2) + 1 ! note int arithmetic rounds TOWARDS zero.
        iend   = count

        oloop : do while (iend > 2)
            if (istart > 1) then ! heap construction
                istart = istart - 1
            else ! heap extraction
                call index_swap(sorted_values(iend),sorted_values(1), sorted_index(iend), sorted_index(1))
                iend = iend - 1
            endif

            iroot = istart
            iloop : do while(iLeftChild(iroot) < iend)
                ichild = iLeftChild(iroot)
                !   Check if there is a right child, then check if that child is greater
                if (ichild < iend .and. sorted_values(ichild) < sorted_values(ichild + 1)) then
                    ichild = ichild + 1
                end if

                if (sorted_values(iroot) < sorted_values(ichild)) then
                    call index_swap(sorted_values(iroot), sorted_values(ichild), sorted_index(iroot), sorted_index(ichild))
                    iroot = ichild
                else
                    exit iloop
                endif

            end do iloop

        end do oloop



    end subroutine heap_sort_index

    pure function iLeftChild(i) result(li)

        implicit none

        integer, intent(in) :: i
        integer :: li

        li = 2 * i

    end function iLeftChild

    pure function iRightChild(i) result(ri)

        implicit none

        integer, intent(in) :: i
        integer :: ri

        ri = (2 * i) + 1

    end function iRightChild

    subroutine index_swap(val1,val2, ind1, ind2)

        implicit none

        real(p2), intent(inout) :: val1, val2
        integer , intent(inout) :: ind1, ind2

        real(p2) :: tval
        integer  :: tind

        tval = val2
        val2 = val1
        val1 = tval

        tind = ind2
        ind2 = ind1
        ind1 = tind

    end subroutine index_swap

    subroutine heapify_index(values,index,count)

        implicit none

        real(p2), dimension(:), intent(inout)  :: values
        integer,  dimension(:), intent(inout)  :: index
        integer,                intent(in   )  :: count


    end subroutine
end module sorting