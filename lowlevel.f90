module lowlevel

    ! this subroutine handles some lowlevel functions, mainly to do with pointers, and reallocating variables.
    ! not sure a better name to call it...

    implicit none

    public :: my_alloc_int_ptr

    contains
!********************************************************************************
    ! This subroutine is useful to expand integer arrays.
    !
    !  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1).
    !  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
    !
    ! For example, with the input
    !
    !  n = 7
    !  x = [9,4,2,1]
    !
    ! this subroutine will return
    !
    !  x = [9,4,2,1,0,0,0].
    !
    ! Note: If n=1, this subroutine takes is as an initialization, and returns
    !
    !  x = [0].
    !
    ! Note: So, this subroutine can only expand an interger array. It does not
    !       shrink an array. If n < size(x), then it is considered as an error and
    !       stop. If you want, I believe you can implement it.
    !
    !********************************************************************************
    subroutine my_alloc_int_ptr(x,n)

    implicit none
    
    integer, intent(in) :: n
    
    integer, dimension(:), pointer :: x
    integer, dimension(:), pointer :: temp
    
    integer :: i
    
    ! Error if n is negative.
    
    if (n <= 0) then
        write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
        stop
    endif
    
    ! Shirinking an array is not implemented... Sorry.
    if ( n < size(x) ) then
        write(*,*) "my_alloc_int_ptr received a smaller dimension. Not implemented. Stop."
        stop
    endif
    
    ! If dimension 1, just allocate and return.
    
    if (n==1) then
        if (associated(x)) nullify(x)
        allocate(x(1))
        return
    endif
    
    ! If reallocation (i.e., n > size(x), create a pointer with a target with the requested dimension.
    
    allocate(temp(n))
    temp = 0
    
    ! Copy the existing data: e.g., for n=7 and x=[9,4,2,1] -> we construct temp=[9,4,2,1,0,0,0].
    
    do i = 1, size(x)
        temp(i) = x(i)
    end do
    
    ! Re-assign the pointer: x=[9,4,2,1,0,0,0].
        x => temp
    
    return
    
    end subroutine my_alloc_int_ptr
    !********************************************************************************

end module lowlevel