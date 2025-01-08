module wall_distance

    use common , only : p2

    implicit none

    private

    public cell_wall_distance

    public compute_wall_distance

    real, dimension(:), allocatable :: cell_wall_distance

    integer, parameter :: ix = 1
    integer, parameter :: iy = 2
    integer, parameter :: iz = 3


    ! Private variables
    integer :: nwall_nodes, sqrt_nwall_nodes
    integer :: nleafs
    integer, dimension(:), allocatable :: wall_nodes
    integer :: ninterior_nodes ! this includes nodes on non-wall boundaries
    integer, dimension(:,:), allocatable :: interior_cells 
    real(p2),dimension(:,:), allocatable :: icell_box_dist
    integer, dimension(:),   allocatable :: tmpinterior_cells
    real(p2),dimension(:),   allocatable :: tmpicell_box_dist

    type bounding_box
        integer                                 :: nwnodes
        integer, dimension(:), pointer          :: wnodes
        integer                                 :: longest_dir
        type(bounding_box), pointer             :: root
        type(bounding_box), pointer             :: branch1
        type(bounding_box), pointer             :: branch2
        logical                                 :: isleaf
        real(p2)                                :: xmin, xmax, ymin, ymax, zmin, zmax
    end type bounding_box

    contains

    subroutine compute_wall_distance

        use grid , only : nnodes, x, y, z, bound, nb, bc_type, ncells, cell

        use sorting , only : heap_sort_index
        
        implicit none

        logical, dimension(nnodes) :: is_wall

        
        real(p2), dimension(:), allocatable :: wnx, wny, wnz, wnx_sortx, wny_sorty, wnz_sortz ! coordinates of wall nodes
        integer,  dimension(:), allocatable :: wns

        type(bounding_box), pointer :: root_box
        ! type(bounding_box)          :: testbox
        type(bounding_box), dimension(:), allocatable :: bbox_leafs

        integer :: ib, inode, iface, ibox, icell
        integer :: ni

        integer :: longest_dir
        integer :: split_number
        real(p2) :: xmin, xmax, ymin, ymax, zmin, zmax
        real(p2) :: dx, dy, dz
        real(p2) :: dummy ! needed for doing some float math with integers...
        
        continue

        ! Define the pointer root_box
        allocate(root_box)

        ! Build the array of interior nodes.
        is_wall = .false.
        bloop1 : do ib = 1,nb
            if (trim(bc_type(ib)) /= 'slip_wall' .and. trim(bc_type(ib)) /= 'no_slip_wall') cycle bloop1
            
            do inode = 1,bound(ib)%nbfaces
                do iface = 2,bound(ib)%bfaces(1,inode) + 1
                    ni = bound(ib)%bfaces(iface,inode)
                    is_wall(ni) = .true.
                end do
            end do
        end do bloop1

        nwall_nodes = count(is_wall)

        ! Compute cutoff (sqrt(nwall_faces))
        dummy = real(nwall_nodes,p2)
        dummy = sqrt(dummy) ! sqrt(int) will return an error
        sqrt_nwall_nodes = int(dummy) 
        nleafs = 0

        allocate(wall_nodes(nwall_nodes))
        allocate(wnx(nwall_nodes))
        allocate(wny(nwall_nodes))
        allocate(wnz(nwall_nodes))
        allocate(wns(nwall_nodes))
        ! allocate(wn_sorty(nwall_nodes))
        ! allocate(wn_sortz(nwall_nodes))

        nwall_nodes = 0
        is_wall = .false.
        bloop2 : do ib = 1,nb
            if (trim(bc_type(ib)) /= 'slip_wall' .and. trim(bc_type(ib)) /= 'no_slip_wall') cycle bloop2
            
            do inode = 1,bound(ib)%nbfaces
                add_node_loop : do iface = 2,bound(ib)%bfaces(1,inode) + 1
                    ni = bound(ib)%bfaces(iface,inode)
                    if (is_wall(ni)) cycle add_node_loop
                    is_wall(ni) = .true.
                    nwall_nodes = nwall_nodes + 1
                    wall_nodes(nwall_nodes) = ni
                    wnx(nwall_nodes)        = x(ni)
                    wny(nwall_nodes)        = y(ni)
                    wnz(nwall_nodes)        = z(ni)
                end do add_node_loop
            end do
        end do bloop2

        ! Set the split direction for the first bounding box
        call sort_longest(nwall_nodes,wnx,wny,wnz,wall_nodes,root_box,wns)

        call construct_bounding_box(nwall_nodes,wall_nodes,root_box)

        allocate(bbox_leafs(nleafs))

        nleafs = 0
        call extract_leafs(nleafs, root_box, bbox_leafs)

        deallocate(root_box)

        ! now we compute the distance to each box
        allocate(interior_cells(nleafs, ncells))
        allocate(icell_box_dist(nleafs, ncells))
        allocate(tmpinterior_cells(nleafs))
        allocate(tmpicell_box_dist(nleafs))

        nloop : do icell = 1,ncells
            if (is_wall(icell)) then
                interior_cells(:,icell) = 0
                cycle nloop
            endif
            do ibox = 1,nleafs
                tmpinterior_cells(ibox) = ibox
                tmpicell_box_dist(ibox) = distance_to_block(bbox_leafs(ibox)%xmin, bbox_leafs(ibox)%xmax, &
                                                               bbox_leafs(ibox)%ymin, bbox_leafs(ibox)%ymax, & 
                                                               bbox_leafs(ibox)%zmin, bbox_leafs(ibox)%zmax, &
                                                               cell(icell)%xc, cell(icell)%yc, cell(icell)%zc )
            end do
            ! fortran doesn't like it when you alias variables in subroutines.  I suppose I could rewrite the heap sort algorithm
            ! with an inout var but I don't want to right now...
            call heap_sort_index(tmpicell_box_dist,tmpinterior_cells,nleafs,interior_cells(:,icell)) 
            do ibox = 1,nleafs
                ib = interior_cells(ibox,icell)
                icell_box_dist(ibox,icell) = tmpicell_box_dist(ib)
            end do

        end do nloop

        deallocate(interior_cells, icell_box_dist)
    end subroutine compute_wall_distance

    recursive subroutine construct_bounding_box(nnodes,nodes,bbox)
    
        use grid , only : x, y, z

        use sorting , only : heap_sort_index
        implicit none

        integer, dimension(:), intent(in) :: nodes  ! the global ID of the wall nodes in this box, sorted in the long dir
        integer,               intent(in) :: nnodes ! number of wall nodes in this bounding box
        ! integer,                        intent(in) :: split  ! index of the last of group 1 (equals nnodes/2)

        ! real(p2), dimension(:),         intent(in) :: wnx,wny,wnz ! x,y, & z coords of the wall nodes
        ! integer,  dimension(:),         intent(in) :: wn_sortx,wn_sorty,wn_sortz

        type(bounding_box), pointer, intent(inout) :: bbox
        type(bounding_box) :: tempbox

        integer :: longest_dir
        integer :: split, nnms
        real(p2) :: xmin, xmax, ymin, ymax, zmin, zmax
        real(p2) :: dx, dy, dz
        integer :: inode, jnode

        ! Branch arrays
        real(p2), dimension(:), allocatable :: wnx1,  wnx2,  wny1,  wny2,  wnz1,  wnz2
        integer , dimension(:), pointer     :: wn1, wn2, wns1, wns2

        bbox%nwnodes = nnodes
        allocate(bbox%wnodes(nnodes))
        bbox%wnodes  = nodes

        nullify(bbox%branch1)
        nullify(bbox%branch2)

        if (bbox%isleaf) then
            nleafs = nleafs + 1
            tempbox = bbox
            ! we don't actually need to do anything else
            return
        end if

        split = nnodes/2
        nnms  = nnodes - split
        allocate(wn1(split), wnx1(split), wny1(split), wnz1(split), wns1(split))
        allocate(wn2(nnms ), wnx2(nnms ), wny2(nnms ), wnz2(nnms ), wns2(nnms ))
        do inode = 1,split
            wn1( inode) =   nodes(inode)
            wnx1(inode) = x(nodes(inode))
            wny1(inode) = y(nodes(inode))
            wnz1(inode) = z(nodes(inode))
        end do
        jnode = 1
        do inode = split+1,nnodes
            wn2(jnode) = nodes(inode)
            wnx2(jnode) = x(nodes(inode))
            wny2(jnode) = y(nodes(inode))
            wnz2(jnode) = z(nodes(inode))
            jnode = jnode + 1
        end do


        ! Repeat recursively for branch 1
        allocate(bbox%branch1)
        bbox%branch1%root => bbox
        call sort_longest(split,wnx1,wny1,wnz1,wn1,bbox%branch1,wns1)

        ! Create a new bounding box with the new dimensions
        call construct_bounding_box(split, wns1, bbox%branch1)

        ! Repeat recursively for branch 2
        allocate(bbox%branch2)
        bbox%branch2%root => bbox
        call sort_longest(nnms ,wnx2,wny2,wnz2,wn2,bbox%branch2,wns2)

        ! Create a new bounding box with the new dimensions
        bbox%branch2%root = bbox
        call construct_bounding_box(nnms, wns2, bbox%branch2)

        deallocate(wn1, wnx1, wny1, wnz1)
        deallocate(wn2, wnx2, wny2, wnz2)
        
    end subroutine construct_bounding_box

    recursive subroutine extract_leafs(nlf, root_box, leaf_array)

        implicit none

        type(bounding_box),               intent(inout) :: root_box
        integer,                          intent(inout) :: nlf
        type(bounding_box), dimension(:), intent(inout) :: leaf_array

        if (root_box%isleaf) then
            nlf = nlf + 1
            leaf_array(nlf) = root_box
            if (associated(leaf_array(nlf)%root)) nullify(leaf_array(nlf)%root)
        else
            call extract_leafs(nlf, root_box%branch1, leaf_array)
            call extract_leafs(nlf, root_box%branch2, leaf_array)
        end if

        if (associated(root_box%root)) nullify(root_box%root)
        if (associated(root_box%branch1)) nullify(root_box%branch1)
        if (associated(root_box%branch2)) nullify(root_box%branch2)
        if (associated(root_box%wnodes))  nullify(root_box%wnodes)

    end subroutine extract_leafs


    subroutine sort_longest(n, wnx, wny, wnz, wn, box, wns)

        use sorting , only : heap_sort_index

        implicit none 

        integer,                intent(in) :: n
        real(p2), dimension(:), intent(in) :: wnx, wny, wnz
        integer,  dimension(:), intent(in) :: wn

        type(bounding_box),     intent(out):: box
        integer,  dimension(:), intent(out):: wns

        real(p2) :: dx, dy, dz

        ! calculate the domain of the wall boundaries for branch 1
        box%xmin = minval(wnx)
        box%xmax = maxval(wnx)
        box%ymin = minval(wny)
        box%ymax = maxval(wny)
        box%zmin = minval(wny)
        box%zmax = maxval(wnz)

        if (n <= sqrt_nwall_nodes) then
            ! no more splitting needed.
            box%isleaf = .true.
            ! we don't need to sort them in this case so instead we just send them back as is.
            wns = wn
            return
        end if

        box%isleaf = .false.

        dx = box%xmax - box%xmin
        dy = box%ymax - box%ymin
        dz = box%zmax - box%zmin

        box%longest_dir = compute_longest_direction(dx,dy,dz)

        ! sort the nodes in the longest direction
        select case(box%longest_dir)
        case(ix)
            call heap_sort_index(wnx, wn, n, wns)
        case(iy)
            call heap_sort_index(wny, wn, n, wns)
        case(iz)
            call heap_sort_index(wnz, wn, n, wns)
        end select

    end subroutine sort_longest

    subroutine split_long(n,split, wxs,wns, wxs1,wxs2, wns1,wns2, groupID)

        implicit none

        integer,                intent(in) :: n, split ! Length of node array, split point
        real(p2), dimension(:), intent(in) :: wxs ! sorted coordinates of wall nodes
        integer , dimension(:), intent(in) :: wns ! pointer of sorted coordinates => global index

        real(p2), dimension(:), intent(out):: wxs1, wxs2 ! sorted split coordinates
        integer , dimension(:), intent(out):: wns1, wns2 ! split sorted pointers
        integer , dimension(:), intent(out):: groupID    ! identifier of which split group each point belongs to

        wxs1 = wxs(1:split)
        wxs2 = wxs(split + 1:n)

        wns1 = wns(1:split)
        wns2 = wns(split + 1:n)

        groupID(1:split)   = 1
        groupID(split+1:n) = 2

    end subroutine split_long

    pure function compute_longest_direction(dx,dy,dz) result(dir)

        implicit none

        real(p2), intent(in) :: dx, dy, dz
        integer :: dir

        real(p2) :: dmax

        dir = ix
        dmax = dx

        if (dy > dmax) then
            dmax = dy
            dir = iy
        end if

        if (dz > dmax) then
            dir = iz
        end if
        
    end function compute_longest_direction

    pure function distance_to_block(xmin,xmax,ymin,ymax,zmin,zmax,px,py,pz) result(dist)

        ! A rather clever method for quickly finding the distance to the closest point along the surface of an axis aligned box.
        ! "borrowed" from: 
        ! https://stackoverflow.com/questions/5254838/calculating-distance-between-a-point-and-a-rectangular-box-nearest-point

        use common , only : zero

        implicit none

        real(p2), intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax ! dimensions of axis aligned box
        real(p2), intent(in) :: px,py,pz ! location of point

        real(p2) :: dist

        real(p2) :: dx,dy,dz

        dx = max(xmin - px,zero,px - xmax)
        dy = max(ymin - py,zero,py - ymax)
        dz = max(zmin - pz,zero,pz - zmax)

        dist = sqrt(dx**2 + dy**2 + dz**2)
    end function distance_to_block
end module wall_distance