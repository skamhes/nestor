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

        integer, dimension(:), allocatable :: base_cell, working_cell, next_cell
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
            do jcell = 1,bound(ib)%nbfaces
                cb = bound(ib)%bcell(jcell)
                if (li_cell(cb) == 0) then
                    nlines = nlines + 1
                    li_cell(cb) = -nlines
                    li_nodes(bound(ib)%bfaces(2:,jcell)) = .true.
                end if
            end do
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
            lines(jline)%ncells = 0
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
            lines(jline)%ncells = 0
        end do


    end subroutine build_lines

end module limplicit