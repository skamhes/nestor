module grid
    use common , only : p2

    implicit none

    public
    ! GRID DATA
    !------------------------------------------
    !>> Node data
    integer                             :: nnodes
    real(p2), dimension(:  )  , pointer :: x, y, z

    !------------------------------------------
    !>> Boundary element data
    integer                              :: nb      !# of boundary segments
    !integer      , dimension(:), pointer :: nbnodes !# of boundary nodes
    !integer      , dimension(:), pointer ::  bnode  !List of boundary nodes
    character(80), dimension(:), pointer :: bc_type !type of boundary condition

    !------------------------------------------
    !>> Element conenctivity data

    !>> Triangular element data
    integer                           :: ntria !# of triangles
    integer , dimension(:,:), pointer ::  tria !List of vertices

    !>> Quadrilateral element data
    integer                           :: nquad !# of quadrilaterals
    integer , dimension(:,:), pointer ::  quad !List of vertices
    
    !>> Tetrahedral element data
    integer                           ::  ntet !# of quadrilaterals
    integer , dimension(:,:), pointer ::   tet !List of vertices

    !>> Tetrahedral element data
    integer                           ::  npyr !# of quadrilaterals
    integer , dimension(:,:), pointer ::   pyr !List of vertices

    !>> Tetrahedral element data
    integer                           ::  nprs !# of quadrilaterals
    integer , dimension(:,:), pointer ::   prs !List of vertices

    !>> Tetrahedral element data
    integer                           ::  nhex !# of quadrilaterals
    integer , dimension(:,:), pointer ::   hex !List of vertices

    contains

    subroutine read_grid

        !********************************************************************************
        !* Read the grid and boundary condition file.
        !* Note: This is a combo of the 2D cell centered and 3D node centered grids
        !        which will be fun...
        !* Note: Currently, the EDU3D-Euler works only for pure tetrahedral grids.
        !*       I think you can modify the code to make it work for other elements.
        !*
        !********************************************************************************
        !*
        !* 1. "datafile_grid_in": .ugrid file name. 
        !*
        !*
        !* 2. "datafile_bcmap_in" is assumed have been written in the following format:
        !*
        !*   -----------------------------------------------------------------------
        !*    write(*,*) nb (# of boundaries)
        !*    do i = 1, nb
        !*     write(*,*) i, bc_type
        !*    end do
        !*   -----------------------------------------------------------------------
        !*
        !*   NOTE: bc_type is the name of the boundary condition, e.g.,
        !*
        !*         1. "freestream"
        !*             Roe flux with freestream condition on the right state.
        !*
        !*         2. "outflow_supersonic"
        !*             Roe flux with the interior state as the right state.
        !*             (equivalent to the interior-extrapolation condition.)
        !*
        !*         3. "slip_wall"
        !*             Right state with zero mass flux through the boundary.
        !*
        !*         4. "outflow_subsonic_p0"
        !*             Fix the back pressure. This should work for subsonic flows in a
        !*             large enough domain.
        !*
        !*         !You can implement more BCs.
        !*
        !*
        !********************************************************************************
        !* Data to be read and stored:
        !*
        !* 1. Some numbers
        !*    nnodes        = number of nodes
        !*    ntria         = number of triangular boundary elements
        !*    nquad         = number of quadrilateral boundary elements
        !*    ntet          = number of tetrahedra
        !*    npyr          = number of pyramids
        !*    nprs          = number of prims
        !*    nhex          = number of hexahedra
        !*
        !* 2. Boundary element data:
        !*    tria = list of vertices of each triangular boundary element
        !*    quad = list of vertices of each quadrilateral boundary element
        !*
        !* 3. Volume element data:
        !*    tet  = list of vertices of each tetrahedron
        !*    pyr  = list of vertices of each pyramid
        !*    prs  = list of vertices of each prism
        !*    hex  = list of vertices of each hexehedron
        !*
        !* 4. x, y, z coordinates of ndoes
        !*    x    = x-coordinate of the nodes
        !*    y    = y-coordinate of the nodes
        !*    z    = z-coordinate of the nodes
        !*
        !* 5. Boundary Data:
        !*    nb      = number of boundary groups
        !*    bc_type = boundary condition name
        !*
        !********************************************************************************
        use inout , only : filename_grid, filename_bc
        implicit none

        integer :: os
        integer :: i, ncells, dummy_int
        !integer , dimension(100,8) ::   dummy_debug ! use to debug public variables
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading : ", trim(filename_grid)
        write(*,*)

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! 1. Read the grid file.

        ! Open the input file.
        open(unit=1, file=filename_grid, status="old", iostat=os)

        if (os .NE. 0) then
            write(*,*)  "FILE ERROR: STOP!"
            stop
        end if
        ! Read: get the size of the grid
        read(1,*) nnodes, ntria, nquad, ntet, npyr, nprs, nhex

        ! Write out the grid data.

        write(*,*) " Total grid numbers:"
        write(*,*) "      Nodes = ", nnodes
        write(*,*) "  Triangles = ", ntria
        write(*,*) "      Quads = ", nquad
        write(*,*) "      Tetra = ", ntet
        write(*,*) "      Hexa  = ", nhex
        write(*,*) "   Pyramids = ", npyr
        write(*,*) "     Prisms = ", nprs
        write(*,*)

        ! Allocate node and element arrays

        if (ntria > 0) allocate(tria(ntria,4))
        if (nquad > 0) allocate(quad(nquad,5))
        if (ntet > 0)  allocate(tet( ntet, 4))
        if (npyr > 0)  allocate(pyr( npyr, 5))
        if (nprs > 0)  allocate(prs( nprs, 6))
        if (nhex > 0)  allocate(hex( nhex, 8))
        
        allocate(x(nnodes),y(nnodes),z(nnodes))

        ! READ: Read the nodal coordinates

        write(*,"(A)",advance="no") " Reading nodes..."

        node_read: do i = 1, nnodes
            read(1,*) x(i), y(i), z(i)
        end do node_read
        write(*,"(A)") "...done"

        write(*,"(A)",advance="no") " Building mesh."
        ! Read element-connectivity information
        ! NOTE: EVERY CELL TYPE MUST BE READ IN A SPECIFIC ORDER

        ! Triangles: assumed that the vertices are ordered counterclockwise
        !
        !         v3
        !         /\
        !        /  \
        !       /    \
        !      /      \
        !     /        \
        !    /__________\
        !   v1           v2

        ! READ: read connectivity info for triangles
        ncells = 1
        if ( ntria > 0) then
            do
                if (ncells > ntria) EXIT
                read(1,*) tria(ncells,1),tria(ncells,2),tria(ncells,3)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "." ! Progress bar...

        ! Quads: assumed that the vertices are ordered counterclockwise
        !
        !        v4________v3
        !         /        |
        !        /         |
        !       /          |
        !      /           |
        !     /            |
        !    /_____________|
        !   v1             v2

        ! READ: read connectivity info for quadrilaterals
        ncells = 1
        if ( nquad > 0) then
            do
                if (ncells > nquad) EXIT
                read(1,*) quad(ncells,1),quad(ncells,2),quad(ncells,3),quad(ncells,4)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."
        
        ! Read Tria Boundary Group Number
        ncells = 1
        if ( ntria > 0) then
            do
                if (ncells > ntria) EXIT
                read(1,*) tria(ncells,4)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."
        
        ! Read Quad Boundary Group Number
        ncells = 1
        if ( nquad > 0) then
            do
                if (ncells > nquad) EXIT
                read(1,*) quad(ncells,5)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."

        ! Read connectivity info for tets
        ncells = 1
        if ( ntet > 0) then
            do
                if (ncells > ntet) EXIT
                read(1,*) tet(ncells,1),tet(ncells,2),tet(ncells,3),tet(ncells,4)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."

        ! Read connectivity info for pyrs
        ncells = 1
        if ( npyr > 0) then
            do
                if (ncells > npyr) EXIT
                read(1,*) pyr(ncells,1),pyr(ncells,2),pyr(ncells,3),pyr(ncells,4), &
                pyr(ncells,5)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."

        ! Read connectivity info for prisms
        ncells = 1
        if ( nprs > 0) then
            do
                if (ncells > nprs) EXIT
                read(1,*) prs(ncells,1),prs(ncells,2),prs(ncells,3),prs(ncells,4), &
                prs(ncells,5),prs(ncells,6)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)",advance="no") "."            

        ! Read connectivity info for hexa
        ncells = 1
        if ( nhex > 0) then
            do
                if (ncells > nhex) EXIT
                read(1,*) hex(ncells,1),hex(ncells,2),hex(ncells,3),hex(ncells,4), &
                hex(ncells,5),hex(ncells,6),hex(ncells,7),hex(ncells,8)
                ncells = ncells + 1
            end do
        end if
        write(*,"(A)") "..done"

        ! Close the file
        close(1)

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! 2. Read the boundary condition data file
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading the boundary condition file: ", &
            trim(filename_bc)
        write(*,*)

        open(unit=2, file=filename_bc, status="unknown", iostat=os)
        read(2,*) nb
        
        allocate(bc_type(nb))
        do i = 1, nb
            read(2,*) dummy_int, bc_type(i)                               
        end do
        !  Print the data
        do i = 1, nb
            write(*,'(a10,i3,a12,a20)') " boundary", i, "  bc_type = ", trim(bc_type(i))
        end do

        write(*,*)

        close(2)
        ! End of Read the boundary condition data file
        !--------------------------------------------------------------------------------

        write(*,*)
        write(*,*) " Finished Reading : ", trim(filename_grid), " and ",&
            trim(filename_bc)
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    end subroutine read_grid

    subroutine read_su2
        ! Read .su2 grid format
        ! If anyone from the su2 group happens to read this maybe don't.  I wrote this at 2am and spent the minimum effort in order
        ! to get it working.  You have been warned........ :)
        use inout, only : filename_grid, filename_bc
        implicit none
        
        integer :: os
        integer :: i, ncells, dummy_int, comment_ind
        integer :: ndim, mark_counter
        integer, dimension(:), allocatable   :: cell_type
        integer, dimension(:,:), allocatable :: input_cells
        character(80) :: buffer, rawbuffer
        logical :: read_dim = .false.
        logical :: read_elm = .false.
        logical :: read_bnd = .false.
        logical :: read_nds = .false.

        type temp_bound
            integer                             ::    nb_faces, nb_tria, nb_quad
            integer, dimension(:,:),allocatable ::    btria 
            integer, dimension(:,:),allocatable ::    bquad 
        end type temp_bound

        type(temp_bound), dimension(:), allocatable :: btemp

        
        !integer , dimension(100,8) ::   dummy_debug ! use to debug public variables
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading : ", trim(filename_grid)
        write(*,*)

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! 1. Read the grid file.

        ! Open the input file.
        open(unit=1, file=filename_grid, status="unknown", iostat=os)
          
        
        nnodes = 0
        ntria = 0
        nquad = 0
        ntet = 0
        npyr = 0
        nprs = 0
        nhex = 0

        do
            read(1,'(A)') rawbuffer
            rawbuffer = adjustl(rawbuffer)
            ! Remove leading white spaces
            comment_ind = index(rawbuffer,'%')
            if (comment_ind == 0) then ! no comment present
                buffer = rawbuffer
            elseif (comment_ind == 1) then ! comment as first char.  Cycle
                cycle
            else
                buffer = rawbuffer(1:comment_ind)
            endif

            if ( index(buffer,'NDIME') .ne. 0 ) then
                backspace(1) ! re-read the line
                read(1,*) rawbuffer, ndim
                if (ndim == 3) then
                    read_dim = .true.
                    cycle ! correct number of dimensions
                elseif (ndim == 2) then
                    write(*,*) '2D flow not supported. Stop!'
                    stop
                else
                    write(*,*) 'Invalid dimension. Stop!'
                    stop
                endif
                cycle
            endif

            if ( index(buffer,'NELEM') .ne. 0 ) then
                backspace(1) ! re-read the line
                read(1,*) rawbuffer, ncells
                allocate(cell_type(ncells))
                allocate(input_cells(8,ncells))
                count_cells : do i = 1,ncells
                    read(1,*) cell_type(i),input_cells(:,i)
                    select case (cell_type(i))
                        case (3) ! line
                            write(*,*) "There shouldn't be a line here.  Something is wrong. Stop!"
                            stop
                        case (5) ! triangle
                            write(*,*) "There shouldn't be a triangle here.  Something is wrong. Stop!"
                            stop
                        case (9)
                            write(*,*) "There shouldn't be a quad here.  Something is wrong. Stop!"
                            stop
                        case (10)
                            ntet = ntet + 1
                        case (12)
                            nhex = nhex + 1
                        case (13)
                            nprs = nprs + 1
                        case (14)
                            npyr = npyr + 1
                        case default
                            write(*,*) "Invalid cell type.  Cell: ", i, ".  Something is wrong. Stop!"
                            stop
                    end select
                end do count_cells
                
                input_cells = input_cells + 1 ! su2 starts counting at 0, I want to start with 1. This kinda feels like a hack...
                ! I just realized all the cells are in row major.  But fixing it is gonna be a pita.  And this is only for loading.
                ! So i guess I'll do it later... famous last words...
                if (ntet > 0)  allocate(tet( ntet, 4))
                if (npyr > 0)  allocate(pyr( npyr, 5))
                if (nprs > 0)  allocate(prs( nprs, 6))
                if (nhex > 0)  allocate(hex( nhex, 8))
                ntet = 0
                npyr = 0
                nprs = 0
                nhex = 0
                add_cells : do i = 1,ncells
                    select case (cell_type(i))
                        case (10)
                            ntet = ntet + 1
                            tet(ntet,:) = input_cells(1:4,i)
                        case (12)
                            nhex = nhex + 1
                            hex(nhex,:) = input_cells(1:8,i)
                        case (13)
                            nprs = nprs + 1
                            prs(nprs,:) = input_cells(1:6,i)
                        case (14)
                            npyr = npyr + 1
                            pyr(npyr,:) = input_cells(1:5,i)
                        case default
                            write(*,*) "Invalid cell type.  Cell: ", i, ".  Something is wrong. Stop!"
                            stop
                    end select
                end do add_cells
                deallocate(cell_type)
                deallocate(input_cells)
                read_elm = .true.
                cycle
            end if

            if ( index(buffer,'NPOIN') .ne. 0 ) then
                backspace(1) ! re-read the line
                read(1,*) rawbuffer, nnodes
                allocate(x(nnodes),y(nnodes),z(nnodes))
                read_nodes : do i = 1,nnodes
                    read(1,*) x(i), y(i), z(i)
                end do read_nodes
                read_nds = .true.
                cycle
            end if

            if ( index(buffer,'NMARK') .ne. 0 ) then
                backspace(1) ! re-read the line
                read(1,*) rawbuffer, nb
                mark_counter = 1
                allocate(btemp(nb))
                loop_marks : do
                    read(1,'(A)') buffer
                    if ( index(buffer,'MARKER_ELEMS') .ne. 0 ) then
                        backspace(1)
                        read(1,*) rawbuffer, btemp(mark_counter)%nb_faces
                        allocate(cell_type(    btemp(mark_counter)%nb_faces))
                        allocate(input_cells(4,btemp(mark_counter)%nb_faces))
                        btemp(mark_counter)%nb_tria = 0
                        btemp(mark_counter)%nb_quad = 0
                        count_bfaces : do i = 1,btemp(mark_counter)%nb_faces
                            read(1,*) cell_type(i),input_cells(:,i)
                            select case (cell_type(i))
                                case (3) ! line
                                    write(*,*) "There shouldn't be a line here.  Something is wrong. Stop!"
                                    stop
                                case (5) ! triangle
                                    btemp(mark_counter)%nb_tria = btemp(mark_counter)%nb_tria + 1
                                case (9)
                                    btemp(mark_counter)%nb_quad = btemp(mark_counter)%nb_quad + 1
                                case default
                                    write(*,*) "Invalid bface type.  Bface: ", i, ".  Something is wrong. Stop!"
                                    stop
                            end select
                        end do count_bfaces
                        input_cells = input_cells + 1
                        if (btemp(mark_counter)%nb_tria > 0)  allocate(btemp(mark_counter)%btria( btemp(mark_counter)%nb_tria, 4))
                        if (btemp(mark_counter)%nb_quad > 0)  allocate(btemp(mark_counter)%bquad( btemp(mark_counter)%nb_quad, 5))
                        btemp(mark_counter)%nb_tria = 0
                        btemp(mark_counter)%nb_quad = 0
                        add_bfaces : do i = 1,btemp(mark_counter)%nb_faces
                            select case (cell_type(i))
                                case (5)
                                    btemp(mark_counter)%nb_tria = btemp(mark_counter)%nb_tria + 1
                                    btemp(mark_counter)%btria(btemp(mark_counter)%nb_tria,1:3) = input_cells(1:3,i)
                                    btemp(mark_counter)%btria(btemp(mark_counter)%nb_tria,  4) = mark_counter
                                case (9)
                                    btemp(mark_counter)%nb_quad = btemp(mark_counter)%nb_quad + 1
                                    btemp(mark_counter)%bquad(btemp(mark_counter)%nb_quad,1:4) = input_cells(1:4,i)
                                    btemp(mark_counter)%bquad(btemp(mark_counter)%nb_quad,5)   = mark_counter
                                case default
                                    write(*,*) "Invalid bface type.  Bface: ", i, ".  Something is wrong. Stop!"
                                    stop
                            end select
                        end do add_bfaces
                        deallocate(input_cells,cell_type)
                    else
                        cycle loop_marks
                    end if
                    mark_counter = mark_counter + 1
                    if (mark_counter > nb) then
                        exit loop_marks
                    end if
                end do loop_marks

                combine_marks1 : do i = 1,nb
                    ntria = ntria + btemp(i)%nb_tria
                    nquad = nquad + btemp(i)%nb_quad
                end do combine_marks1

                if (ntria > 0) allocate(tria(ntria,4))
                if (nquad > 0) allocate(quad(nquad,5))
                ntria = 0
                nquad = 0
                combine_marks2 : do i = 1,nb
                    if (btemp(i)%nb_tria > 0) then
                        tria(ntria+1:(ntria + btemp(i)%nb_tria),:) = btemp(i)%btria(:,:)
                        deallocate(btemp(i)%btria)
                    end if
                    if (btemp(i)%nb_quad > 0) then
                        quad(nquad+1:nquad + btemp(i)%nb_quad,:) = btemp(i)%bquad(:,:)
                        deallocate(btemp(i)%bquad)
                    end if
                    ntria = ntria + btemp(i)%nb_tria
                    nquad = nquad + btemp(i)%nb_quad
                    
                end do combine_marks2
                deallocate(btemp)
                read_bnd = .true.
            end if

            if (read_bnd .and. read_dim .and. read_elm .and. read_nds) exit ! everything is done
        end do

        close(1)

        write(*,*) ' Finished reading .su2 grid.'
        write(*,*)
        write(*,*) " Total grid numbers:"
        write(*,*) "      Nodes = ", nnodes
        write(*,*) "  Triangles = ", ntria
        write(*,*) "      Quads = ", nquad
        write(*,*) "      Tetra = ", ntet
        write(*,*) "      Hexa  = ", nhex
        write(*,*) "   Pyramids = ", npyr
        write(*,*) "     Prisms = ", nprs
        write(*,*)


         !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! 2. Read the boundary condition data file
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Reading the boundary condition file: ", &
            trim(filename_bc)
        write(*,*)

        open(unit=2, file=filename_bc, status="unknown", iostat=os)
        read(2,*) nb
        
        allocate(bc_type(nb))
        do i = 1, nb
            read(2,*) dummy_int, bc_type(i)                               
        end do
        !  Print the data
        do i = 1, nb
            write(*,'(a10,i3,a12,a20)') " boundary", i, "  bc_type = ", trim(bc_type(i))
        end do

        write(*,*)

        close(2)
        ! End of Read the boundary condition data file
        !--------------------------------------------------------------------------------

        write(*,*)
        write(*,*) " Finished Reading : ", trim(filename_grid), " and ",&
            trim(filename_bc)
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    end subroutine read_su2
end module grid