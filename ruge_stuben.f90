module ruge_stuben


    ! use common , only : p2

    implicit none

    private
    
    integer, parameter :: RS_U = 0 ! unassigned
    integer, parameter :: RS_C = 1
    integer, parameter :: RS_F = 2
    integer, parameter :: RS_Fnew = 3

    ! logical, parameter :: STRONG_AGGLOM = .true.
    
    ! real(p2), parameter    :: THRESHOLD = 0.25_p2
    
    public rs_agglom
    
    contains

    subroutine rs_agglom(ncells,level,C,R,ngroups,restrictR,restrictC,prolongR,prolongC)

        implicit none

        integer, intent(in)                     :: ncells, level
        ! real(p2),dimension(:,:,:), intent(in)   :: V
        integer, dimension(:), intent(in)       :: C
        integer, dimension(:), intent(in)       :: R

        integer, dimension(:), pointer, intent(out)  :: restrictR
        integer, dimension(:),          intent(out)  :: restrictC
        integer, dimension(:),          intent(out)  :: prolongR
        integer, dimension(:),          intent(out)  :: prolongC
        integer,                        intent(out)  :: ngroups

        integer, dimension(ncells)      :: w, CF
        integer, dimension(ncells)      :: sorted_to_w  ! vector of weight that corresponds to ith sorted index
        integer                         :: nC
        logical                         :: strong_agglom

        strong_agglom = (level == 1) ! does this work?

        call rs_weight(C,R,ncells, w)

        call rs_build_CF(ncells,C,R,w,CF,nc,sorted_to_w)

        call rs_build_r_p1(ncells,strong_agglom,sorted_to_w,nc,CF,C,R,ngroups,restrictC,restrictR,prolongR,prolongC)

    end subroutine rs_agglom

    subroutine rs_weight(C,R, ncells, w) 
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

        !real(p2),dimension(:,:,:), intent(in)   :: V        ! Value array of jacobian
        integer, dimension(:), intent(in)       :: C        ! column indices of association matrix (non-zero jacobian blocks)
        integer, dimension(:), intent(in)       :: R        ! pointer of row starts in CSR format
        integer              , intent(in)       :: ncells

        ! integer, dimension(:), intent(out)      :: Sr       ! column indices of strongly dependent cells
        ! integer, dimension(:), intent(out)      :: Sc       ! pointer of row starts in CSR format
        integer, dimension(:)                   :: w        ! vector of the number of cells strongly influenced by i

        integer :: i, j
        integer :: cj
        !integer :: nstrong
        !real(p2) :: ti ! threshold_i = |A_ik|*theta

        w(:) = 0
        ! nstrong = 0

        do i = 1,ncells
            ! ti = THRESHOLD * maxval( abs(V(1, 1, R(i):R(i+1)-1)) )
            f_loop : do j = R(i),R(i+1)-1
                cj = c(j)
                w(cj) = w(cj) + 1 ! all face neighbors are automatically added
            end do f_loop 
        end do

    end subroutine rs_weight

    subroutine rs_build_CF(ncells,C,R,w,CF,nC,sorted_to_w)

        implicit none

        integer, dimension(:), target, intent(in)   :: C    ! columns of strongly 
        integer, dimension(:),         intent(in)   :: R
        integer,                       intent(in)   :: ncells
        integer, dimension(:),         intent(inout):: w

        integer, dimension(:),         intent(out)  :: CF ! I briefly considered doing this as a logical but a brief google search 
        ! seems to say that logicals occupy the same amount of space in memory as integers so this is easier to read.
        integer,                       intent(out)  :: nC ! number of C-cells (and number of groups since each group has one C-cell)
        integer, dimension(:)                       :: sorted_to_w  ! vector of weight that corresponds to ith sorted index

        integer                             :: i,j,k
        integer                             :: cj,ck
        integer                             :: wmax
        
        integer, dimension(:), pointer      :: si, sj, sintr
        integer                             :: nintr

        
        integer, dimension(ncells)          :: w_to_sorted  ! vector of sorted index that corresponds to ith
        integer, dimension(:), allocatable  :: wcounter 
        integer, dimension(:), allocatable  :: psorted ! Pointer to first sorted value (same method as CSR)

        integer :: nweights
        integer :: cur_pos, end_pos, tmp, endwi
        integer :: icell
        integer :: wi, vi

        wmax = maxval(w)
        nweights = max(ncells,2*wmax)
        CF = 0
        nc = 0

        allocate(wcounter(0:nweights)) ! Zero index!!! becuse theoretically we could have a weight of 0
        wcounter = 0
        allocate(psorted( 0:nweights+1))
        do icell = 1,ncells
            wi = w(icell)
            wcounter(wi) = wcounter(wi) + 1
        end do

        psorted(0) = 1
        do i = 0,nweights
            psorted(i+1) = psorted(i) + wcounter(i)
            wcounter(i)  = 0 
        end do

        do icell = 1,ncells
            wi = w(icell)
            vi = psorted(wi) + wcounter(wi)
            sorted_to_w(vi)    = icell
            w_to_sorted(icell) = vi
            wcounter(wi) = wcounter(wi) + 1
        end do

        max_loop : do i = nweights,1,-1
            ! Pick out the cell with the max weight
            wi = sorted_to_w(i)
            if (CF(wi) /= 0) cycle max_loop ! already assigned
            CF(wi) = RS_C
            nc = nc + 1
            nloop1 : do j = R(wi),R(wi+1)-1
                cj = C(j)
                if (CF(cj) /= 0) cycle nloop1
                CF(cj) = RS_Fnew
            end do nloop1
            ! Increment the weight on the surrounding cells
            nloopj : do j = R(wi),R(wi+1)-1
                cj = C(j)
                if (CF(cj) /= RS_Fnew) cycle nloopj
                CF(cj) = RS_F
                nloopk : do k = R(cj),R(cj+1)-1
                    ck = C(k)
                    if (CF(ck) /= 0) cycle nloopk
                    w(ck) = w(ck) + 1
                    
                    ! Swap places with the last variable w/ w(ck)'s old weight
                    cur_pos = w_to_sorted(ck) ! current sort index of ck
                    psorted(w(ck)) = min(i-1,psorted(w(ck)) - 1) ! decrease the start index by one
                    end_pos = psorted(w(ck))   ! sort index that ck is moving to, must be less than i to be visited
                    endwi = sorted_to_w(end_pos) ! weight index that is moving to cur_pos
                    ! swap w pointers
                    w_to_sorted(ck)    = end_pos
                    w_to_sorted(endwi) = cur_pos
                    ! swap sorting pointers
                    tmp = sorted_to_w(cur_pos)
                    sorted_to_w(cur_pos) = sorted_to_w(end_pos)
                    sorted_to_w(end_pos) = tmp
                end do nloopk
            end do nloopj
        end do max_loop

        ! ! Ensure every cell is assigned
        ! where (CF == 0) CF = 2
        ! This is not needed as we visit every cell

        ! 2nd Pass
        sp_rs : do i = 1,ncells
            if (CF(i) == RS_C) cycle sp_rs
            si => C(R(i):R(i+1)-1)
            si_loop : do j = 1,(R(i+1)-R(i))
                cj = si(j)
                if ( CF(cj) == RS_C ) cycle si_loop
                sj => C(R(cj):R(cj+1)-1)
                call set_intersect(si,sj,sintr,nintr)
                if (nintr == 0) cycle si_loop
                do k = 1,nintr
                    if (CF(sintr(k)) == RS_C) cycle si_loop
                end do
                CF(cj) = RS_C
                nC = nC + 1
                if (associated(sintr)) deallocate(sintr)
            end do si_loop
        end do sp_rs


    end subroutine rs_build_CF

    subroutine rs_build_r_p1(ncells,strong_agglom,sorted_to_w,nc,CF,C,R,ngroups,restrictC,restrictR,prolongR,prolongC)

        implicit none

        integer,               intent(in)           :: ncells, nc
        logical,               intent(in)           :: strong_agglom
        integer, dimension(:), intent(in)           :: sorted_to_w
        integer, dimension(:), intent(inout)        :: CF ! we will set CF back to unassigned to signal it's been added to P & R   
        integer, dimension(:), intent(in)           :: C
        integer, dimension(:), intent(in)           :: R

        integer,                        intent(out) :: ngroups
        integer, dimension(:), pointer, intent(out) :: restrictR
        integer, dimension(:),          intent(out) :: prolongR
        integer, dimension(:),          intent(out) :: prolongC
        integer, dimension(:),          intent(out) :: restrictC

        integer, dimension(nc+1) :: tmprestrictR
        integer :: i, j, k, l
        integer :: ci, cj, ck, cl
        integer :: ccounter
        
        integer :: strong_agglom_int

        ! initialize some values
        strong_agglom_int = 0
        if (strong_agglom) strong_agglom_int = 1
        ccounter           = 1
        prolongR(:)  = (/ (i, i=1,ncells+1)/) ! 1 value per row
        tmprestrictR(1) = 1
        ngroups      = 0

        cloop : do i = 1,ncells
            ! For now we are starting with the lowest weights.  The idea being they're less likely to "hog" F-cells from other
            ! C-cells
            ci = sorted_to_w(i)
            if ( CF(ci) /= RS_C ) cycle cloop
            ngroups = ngroups + 1
            floop1 : do j = R(ci),R(ci+1)-1
                cj = C(j)
                if (CF(cj) == RS_U) then ! already added
                    cycle floop1
                elseif (CF(cj) == RS_F) then
                    restrictC(ccounter) = cj
                    prolongC( ccounter) = ngroups
                    CF(cj) = RS_Fnew * strong_agglom_int
                    ccounter = ccounter + 1
                elseif (cj == ci) then
                    restrictC(ccounter) = cj
                    prolongC( ccounter) = ngroups
                    CF(cj) = RS_U
                    ccounter = ccounter + 1
                end if
            end do floop1
            if (.not.STRONG_AGGLOM) then 
                tmprestrictR(ngroups+1) = ccounter
                cycle cloop
            endif
            ! With agressive (strong) coarsening we add 2nd face neighbors
            floop2 : do j = R(ci),R(ci+1)-1
                cj = C(j)
                if (CF(cj) /= RS_Fnew) cycle floop2
                CF(cj) = RS_U
                f2loop : do k = R(cj),R(cj+1)-1 ! loop 2nd face neighbors
                    ck = C(k)
                    ! If a C cell is a 2nd face nghbr to ci we will add it to the agglomeration
                    if ( CF(ck) /= RS_C) cycle f2loop
                    floop3 : do l = R(ck),R(ck+1)-1
                        cl = C(l)
                        if (CF(cl) == RS_F) then
                            restrictC(ccounter) = cl
                            prolongC( ccounter) = ngroups
                            CF(cl) = RS_U
                            ccounter = ccounter + 1
                        elseif (cl == ck) then
                            restrictC(ccounter) = cl
                            prolongC( ccounter) = ngroups
                            CF(cl) = RS_U
                            ccounter = ccounter + 1
                        end if
                    end do floop3
                end do f2loop
            end do floop2
        tmprestrictR(ngroups+1) = ccounter
        end do cloop

        allocate(restrictR(ngroups+1))
        restrictR(:) = tmprestrictR(1:ngroups+1)


    end subroutine rs_build_r_p1 

    subroutine rs_build_r_p2(ncells,strong_agglom,sorted_to_w,nc,CF,C,R,ngroups,restrictC,restrictR,prolongR,prolongC)

        ! In this aggressive agglomeration routine nearby C-points need to paths to be included.  So this is a p=2,l=2 method.
        ! Initially CF of all C-points is set to RS_C (1). To test if a C-point has 2 connection paths, the first time it is visited
        ! its CF is set to -i.  On the second pass we check if CF = -i.  If it does we know p=>2 and we can add the point.
        ! As a result, the CF value of C-points can either be RS_C (1), or negative (lt 0).  Is this the best way to do it? Probably
        ! not but I've discovered that interpolation techniques are a massive rabbit hole, and I don't really want to dive down
        ! it right now...

        implicit none

        integer,               intent(in)           :: ncells, nc
        logical,               intent(in)           :: strong_agglom
        integer, dimension(:), intent(in)           :: sorted_to_w
        integer, dimension(:), intent(inout)        :: CF ! we will set CF back to unassigned to signal it's been added to P & R   
        integer, dimension(:), intent(in)           :: C
        integer, dimension(:), intent(in)           :: R

        integer,                        intent(out) :: ngroups
        integer, dimension(:), pointer, intent(out) :: restrictR
        integer, dimension(:),          intent(out) :: prolongR
        integer, dimension(:),          intent(out) :: prolongC
        integer, dimension(:),          intent(out) :: restrictC

        integer, dimension(nc+1) :: tmprestrictR
        integer :: i, j, k, l
        integer :: ci, cj, ck, cl
        integer :: ccounter
        
        integer :: strong_agglom_int

        ! initialize some values
        strong_agglom_int = 0
        if (strong_agglom) strong_agglom_int = 1
        ccounter           = 1
        prolongR(:)  = (/ (i, i=1,ncells+1)/) ! 1 value per row
        tmprestrictR(1) = 1
        ngroups      = 0

        cloop : do i = 1,ncells
            ! For now we are starting with the lowest weights.  The idea being they're less likely to "hog" F-cells from other
            ! C-cells
            ci = sorted_to_w(i)
            if ( CF(ci) /= RS_C .AND. CF(ci) >=0 ) cycle cloop
            ngroups = ngroups + 1
            floop1 : do j = R(ci),R(ci+1)-1
                cj = C(j)
                if (CF(cj) == RS_U) then ! already added
                    cycle floop1
                elseif (CF(cj) == RS_F) then
                    restrictC(ccounter) = cj
                    prolongC( ccounter) = ngroups
                    CF(cj) = RS_Fnew * strong_agglom_int
                    ccounter = ccounter + 1
                elseif (cj == ci) then
                    restrictC(ccounter) = cj
                    prolongC( ccounter) = ngroups
                    CF(cj) = RS_U
                    ccounter = ccounter + 1
                end if
            end do floop1
            if (.not.STRONG_AGGLOM) then 
                tmprestrictR(ngroups+1) = ccounter
                cycle cloop
            endif
            ! With agressive (strong) coarsening we add 2nd face neighbors
            floop2 : do j = R(ci),R(ci+1)-1
                cj = C(j)
                if (CF(cj) /= RS_Fnew) cycle floop2
                CF(cj) = RS_U
                f2loop : do k = R(cj),R(cj+1)-1 ! loop 2nd face neighbors
                    ck = C(k)
                    ! If a C cell is a 2nd face nghbr to ci we will add it to the agglomeration
                    if ( CF(ck) /= RS_C .AND. CF(ci) >=0 ) cycle f2loop
                    floop3 : do l = R(ck),R(ck+1)-1
                        cl = C(l)
                        if (CF(cl) == RS_F) then
                            restrictC(ccounter) = cl
                            prolongC( ccounter) = ngroups
                            CF(cl) = RS_U
                            ccounter = ccounter + 1
                        elseif (cl == ck) then
                            restrictC(ccounter) = cl
                            prolongC( ccounter) = ngroups
                            CF(cl) = RS_U
                            ccounter = ccounter + 1
                        end if
                    end do floop3
                end do f2loop
            end do floop2
        tmprestrictR(ngroups+1) = ccounter
        end do cloop

        allocate(restrictR(ngroups+1))
        restrictR(:) = tmprestrictR(1:ngroups+1)


    end subroutine rs_build_r_p2 

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