module ruge_stuben

    implicit none

    private
    
    integer, parameter :: RS_C = 1
    integer, parameter :: RS_F = 2
    
    public rs_agglom
    
    contains

    subroutine rs_agglom(ncells,C,R,restrictR,restrictC,prolongR,prolongC)

        implicit none

        integer, intent(in)                 :: ncells
        integer, dimension(:), intent(in)   :: C
        integer, dimension(:), intent(in)   :: R

        integer, dimension(:), pointer, intent(out)  :: restrictR
        integer, dimension(:),          intent(out)  :: restrictC
        integer, dimension(:),          intent(out)  :: prolongR
        integer, dimension(:),          intent(out)  :: prolongC

        integer, dimension(ncells)  :: w, CF
        integer                     :: ngroups

        w = rs_weight(C,R,ncells)

        call rs_build_CF(ncells,C,R,w,CF,ngroups)

        call rs_build_r(ncells,ngroups,CF,C,R,restrictC,restrictR,prolongR,prolongC)

    end subroutine rs_agglom

    function rs_weight(C,R, ncells) result(w)
        ! Computes the number of strong influences for each cell i.  Note: ordinarily strong influence would be defined as some 
        ! value |a_ij| > theta, where theta is some threashold.  In this case we are setting the threashold to be any nonzero 
        ! index.  As a result we only need the pointer vector R (in CSR format).
        ! If I decide to change this in the future I can add it to this function.
        ! 
        ! Some very quick analysis has this function as O(ncells*nneighbor^2) where nneighbor <= ncells - 1 for all ncells, and 
        ! usually nneighbor << ncells for very large ncells (for a sparse matrix).  This means in the limit as ncells --> inf it is
        ! ~O(ncells).
        !
        ! I think...

        implicit none

        integer, dimension(:)       :: C        ! column indices of strongly dependent cells
        integer, dimension(:)       :: R        ! pointer of row starts in CSR format
        integer                     :: ncells

        integer, dimension(ncells)  :: w        ! vector of the number of cells strongly influenced by i

        integer :: i, j
        integer :: cj

        w(:) = 0

        do i = 1,ncells
            do j = R(i),R(i+1)-1
                cj = c(j)
                w(cj) = w(cj) + 1
            end do
        end do

    end function rs_weight

    subroutine rs_build_CF(ncells,C,R,w,CF,nC)

        implicit none

        integer, dimension(:), target, intent(in)   :: C
        integer, dimension(:),         intent(in)   :: R
        integer,                       intent(in)   :: ncells
        integer, dimension(:),         intent(inout):: w

        integer, dimension(:),         intent(out)  :: CF ! I briefly considered doing this as a logical but a brief google search 
        ! seems to say that logicals occupy the same amount of space in memory as integers so this is easier to read.
        integer,                       intent(out)  :: nC ! number of C-cells (and number of groups since each group has one C-cell)

        integer                             :: nadded
        logical, dimension(ncells)          :: isadded
        integer                             :: i,j,k
        integer                             :: ni, nn
        integer                             :: cj,ck
        integer                             :: wmax
        
        integer, dimension(:), pointer      :: si, sj, sintr
        integer                             :: nintr

        nullify(sintr)

        ! First pass
        isadded = .false.
        nadded  = 0 
        nc = 0

        ! define a start variable for our search
        i = 1

        fp_rs : do while ( nadded < ncells )
            do while ( isadded(i) ) ! breaks when isadded(i) = .false.
                ! This is slightly convoluted.  I will only increment if isadded(i) = .true. which means we should not accidentally
                ! pass an unadded cell.  This means we don't have to start at i=1 each time we cycle the outer loop.
                i = i + 1
            end do

            ! i is now our start vertex
            ! Travel along nodes until we find a local maximum of w.
            wmax = -1
            ni = i
            nn = -1
            do while (ni /= nn)
                nloop1 : do j = R(ni),R(ni+1)-1
                    cj = C(j)
                    if (isadded(cj)) cycle nloop1
                    if (w(cj) > wmax) then
                        wmax = w(cj)
                        nn = cj
                    endif
                end do nloop1
            end do

            ! Now ni represents a local maxima
            ! assign it to {C}
            CF(ni) = RS_C
            nC = nC + 1
            isadded(ni) = .true.
            nadded = nadded + 1
            
            ! Assign all the unassigned neighbors to {F}
            nloop2 : do j = R(ni),R(ni+1)-1
                cj = C(j)
                if (isadded(cj)) cycle nloop2
                CF(cj) = RS_F
                isadded(cj) = .true.
                nadded = nadded + 1
                ! Increment weight of neighbors to nodes in {F_new}
                nk_loop : do k = R(cj),R(cj+1)-1
                    ! technically this could increment the weight of cells that are about to be added, but once they're added their
                    ! weight doesn't matter and creating a seperate vector of neighbors seems more expensive.
                    ck = C(k)
                    if (isadded(ck)) cycle nk_loop
                    w(ck) = w(ck) + 1
                end do nk_loop 
            enddo nloop2
        end do fp_rs

        ! 2nd Pass
        sp_rs : do i = 1,ncells
            if (CF(i) == RS_C) cycle sp_rs
            si => C(R(i):R(i+1)-1)
            si_loop : do j = 1,(R(i+1)-R(i))
                cj = si(i)
                if ( CF(cj) == RS_C ) cycle si_loop
                sj => C(R(cj):R(cj+1)-1)
                call set_intersect(si,sj,sintr,nintr)
                if (nintr == 0) cycle si_loop
                CF(cj) = RS_C
                nC = nC + 1
                if (associated(sintr)) deallocate(sintr)
            end do si_loop
        end do sp_rs

    end subroutine rs_build_CF

    subroutine rs_build_r(ncells,ngroups,CF,C,R,restrictC,restrictR,prolongR,prolongC)

        implicit none

        integer,               intent(in)           :: ncells, ngroups
        integer, dimension(:), intent(in)           :: CF   
        integer, dimension(:), intent(in)           :: C
        integer, dimension(:), intent(in)           :: R

        integer, dimension(:), pointer, intent(out) :: restrictR
        integer, dimension(:),          intent(out) :: prolongR
        integer, dimension(:),          intent(out) :: prolongC
        integer, dimension(:),          intent(out) :: restrictC

        logical, dimension(ncells) :: added
        integer :: igroup, ic, j
        integer :: cj
        integer :: ccounter
        
        ! initialize some values
        ic           = 1
        allocate(restrictR(ngroups))
        prolongR(:) = (/ (ic, ic=1,ncells+1)/) ! 1 value per row
        restrictR(1) = 1
        ccounter     = 1
        added        = .false.

        do igroup = 1,ngroups
            do while(CF(ic) == RS_F)
                ic = ic + 1
            end do

            rloop : do j = R(ic),R(ic+1)-1
                cj = C(j)
                if (added(cj)) cycle rloop
                if (cj == ic) then
                    restrictC(ccounter) = cj
                    prolongC(ccounter)  = igroup
                    added(cj) = .true.
                    ccounter = ccounter + 1
                elseif (CF(cj) == RS_C) then
                    cycle rloop
                else ! CF == RS_F
                    restrictC(ccounter) = cj
                    prolongC(ccounter)  = igroup
                    added(cj) = .true.
                    ccounter = ccounter + 1
                endif
            end do rloop

            restrictR((igroup+1)) = ccounter
            

        end do






    end subroutine rs_build_r   

    subroutine set_intersect(si,sj,sintr,nintr)

        implicit none

        integer, dimension(:), intent(in)           :: si, sj ! Input sets to be intersected

        integer, dimension(:), pointer, intent(out) :: sintr
        integer,                        intent(out) :: nintr

        integer :: i, szi, szj
        integer, dimension(:), allocatable :: tempintr

        nintr = 0

        szi = size(si)
        szj = size(sj)
        allocate( tempintr( max(szi,szj) ) )

        do i = 1,szi
            if (.not. any(si(i) == sj) ) cycle
            nintr = nintr + 1
            tempintr(nintr) = i
        end do

        if (nintr > 0) then
            allocate(sintr(nintr))
            sintr = tempintr(1:nintr)
        end if

        deallocate(tempintr)

    end subroutine set_intersect
end module ruge_stuben