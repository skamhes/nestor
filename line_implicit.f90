module limplicit

    implicit none

    public

    integer :: nlines

    ! structure describing a line of cells for line-implicit solver
    type :: line_type
        integer                             :: ncells           ! number of cells in the line
        integer, dimension(:), allocatable  :: lcells           ! vector of the cells in the line
    end type line_type

    type(line_type), dimension(:), allocatable :: lines
    
    contains

    subroutine build_lines

        use grid , only : ncells, cell, nnodes, nb, bound, bc_type

        use utils , only : ibc_type, BC_VISC_STRONG
        implicit none

        integer, dimension(ncells) :: li_cell ! array of the line index of cells
        logical, dimension(nnodes) :: li_nodes ! array of used nodes

        integer, dimension(:), allocatable :: base_cell, working_cell, next_cell, last_cell
        logical, dimension(:), allocatable :: working

        integer, dimension(6) :: candidate_cells
        integer               :: nccells

        integer :: ib, jcell, cb, jline, kcell
        integer :: nworking, cw, cn
        
        ! initialize the two arrays
        li_cell = 0
        li_nodes = .false.

        nlines = 0

        ! count the number of lines
        do ib = 1,nb
            if (ibc_type(ib) /= BC_VISC_STRONG) cycle
            bc_loop : do jcell = 1,bound(ib)%nbfaces
                cb = bound(ib)%bcell(jcell)
                if (.not.(cell(cb)%nvtx == 8 .or. & ! is not hex
                         (cell(cb)%nvtx == 6) ) ) cycle bc_loop! is not prism
                if (li_cell(cb) == 0) then
                    nlines = nlines + 1
                    li_cell(cb) = -nlines
                    li_nodes(bound(ib)%bfaces(2:,jcell)) = .true.
                end if
            end do bc_loop
        end do

        ! allocate the line struct vector
        allocate(lines(nlines))
        allocate(base_cell(nlines))
        allocate(working_cell(nlines))
        allocate(working(nlines))
        allocate(next_cell(nlines))

        working = .true.

        nlines = 0

        ! Set the base and working cells
        do ib = 1,nb
            if (ibc_type(ib) /= BC_VISC_STRONG) cycle
            do jcell = 1,bound(ib)%nbfaces
                cb = bound(ib)%bcell(jcell)
                if (li_cell(cb) < 0) then
                    nlines               = nlines + 1
                    li_cell(cb)          = nlines
                    base_cell(nlines)    = cb
                    working_cell(nlines) = cb
                end if
            end do
        end do

        nworking = nlines ! lines still growing.

        do jline = 1,nlines
            lines(jline)%ncells = 1 ! the base cell has already been counted
        end do
        
        do while(nworking > 0)
            ! mark all the candidate cells
            line_loop : do jline = 1,nlines
                if (.not. working(jline)) cycle line_loop
                ! Identify the candidate cells
                nccells = 0
                cw = working_cell(jline)
                cloop : do kcell = 1,cell(cw)%nnghbrs
                    cn = cell(cw)%nghbr(kcell)
                    if (li_cell(cn) /= 0) cycle cloop ! check if already assigned to line
                    if (.not.(cell(cn)%nvtx == 8 .or. & ! is not hex
                             (cell(cn)%nvtx == 6) ) ) cycle cloop! is not prism
                    if (any(li_nodes(cell(cn)%vtx))) cycle cloop ! new layer => no marked nodes.
                    ! If we've made it this far it's candidate cell
                    nccells = nccells + 1
                    candidate_cells(nccells) = cn
                end do cloop
                if (nccells == 1) then ! only one candidate.  We use it
                    cn = candidate_cells(1)
                    li_cell(cn) = jline
                    next_cell(jline) = cn
                    lines(jline)%ncells = lines(jline)%ncells + 1
                else ! 0 candidate cells means end of the line. 2+ candidates shouldn't happen.  But either way it represents an 
                     ! edge case we're not going to deal with.
                    nworking = nworking - 1
                    working(jline) = .false.
                end if
            end do line_loop
            line_loop2 : do jline = 1,nlines
                if (.not. working(jline)) cycle line_loop2
                cw = working_cell(jline)
                li_nodes(cell(cw)%vtx) = .true. ! we have to wait until we've finished the next line to mark the nodes
                working_cell(jline) = next_cell(jline)
            end do line_loop2
        end do

        do jline = 1,nlines
            allocate(lines(jline)%lcells(lines(jline)%ncells))
            lines(jline)%lcells(1) = base_cell(jline)
            lines(jline)%ncells = 1
        end do

        nworking = nlines
        working  = .true.
        working_cell = base_cell
        last_cell    = base_cell

        ! now assign them to the line structures.
        do while(nworking > 0)
            ! mark all the candidate cells
            line_loop3 : do jline = 1,nlines
                if (.not. working(jline)) cycle line_loop3
                ! Identify the candidate cells
                nccells = 0
                cw = working_cell(jline)
                cloop2 : do kcell = 1,cell(cw)%nnghbrs
                    cn = cell(cw)%nghbr(kcell)
                    if (li_cell(cn) /= jline) cycle cloop2 ! not in this line
                    if (cn == last_cell(jline)) cycle cloop2! wrong direction
                    ! If we've made it this far it's the next cell
                    next_cell(jline) = cn
                    exit cloop2
                end do cloop2
                if (next_cell(jline) == working_cell(jline)) then ! no next cell was found
                    nworking = nworking - 1
                    working(jline) = .false.
                else
                    lines(jline)%ncells = lines(jline)%ncells + 1
                    lines(jline)%lcells(lines(jline)%ncells) = next_cell(jline)
                    last_cell(jline)    = working_cell(jline)
                    working_cell(jline) = next_cell(jline)
                endif
            end do line_loop3
        end do

    end subroutine build_lines

end module limplicit