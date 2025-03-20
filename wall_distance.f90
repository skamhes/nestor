module wall_distance

    use common , only : p2

    implicit none

    private

    public cell_wall_distance

    public compute_wall_distance

    real(p2), dimension(:), allocatable :: cell_wall_distance

    integer, parameter :: ix = 1
    integer, parameter :: iy = 2
    integer, parameter :: iz = 3

    real(p2), parameter :: EXACT_DISTANCE_THRESHOLD = 10.0_p2


    ! Private variables
    integer :: nwall_nodes, sqrt_nwall_nodes
    integer :: nleafs
    integer, dimension(:),   allocatable :: wall_nodes
    integer, dimension(:),   allocatable :: interior_cells 
    real(p2),dimension(:),   allocatable :: icell_box_dist
    integer, dimension(:),   allocatable :: tmpinterior_cells
    real(p2),dimension(:),   allocatable :: tmpicell_box_dist

    type bounding_box
        integer                                 :: nwnodes
        integer, dimension(:), allocatable      :: wnodes
        integer                                 :: longest_dir
        type(bounding_box), pointer             :: root
        type(bounding_box), pointer             :: branch1
        type(bounding_box), pointer             :: branch2
        logical                                 :: isleaf
        real(p2)                                :: xmin, xmax, ymin, ymax, zmin, zmax
    end type bounding_box

    type wnfaces ! boundary faces attached to wall nodes
        integer, dimension(:), allocatable :: bface
        integer, dimension(:), allocatable :: bound
        integer                            :: nfaces
        integer                            :: ni ! mainly for debugging purposes
    end type wnfaces


    contains

    subroutine compute_wall_distance

        use common , only : my_huge

        use grid , only : nnodes, x, y, z, bound, nb, bc_type, ncells, cell

        use sorting , only : heap_sort_index
        
        implicit none

        logical, dimension(nnodes) :: is_wall

        
        real(p2), dimension(:), allocatable :: wnx, wny, wnz ! coordinates of wall nodes
        integer,  dimension(:), allocatable :: wns
        integer,  dimension(:), allocatable :: gnode_to_wnode ! pointer for translating between global nodes and wall nodes 
        integer,  dimension(:), allocatable :: nf ! number of attached wall faces

        type(bounding_box), pointer :: root_box
        ! type(bounding_box)          :: testbox
        type(bounding_box), dimension(:), allocatable :: bbox_leafs

        type(wnfaces), dimension(:), allocatable :: wn_to_wf

        integer :: ib, inode, iface, ibox, icell
        integer :: ni, bxi, fi, wni
        integer :: closest_node

        real(p2) :: dummy ! needed for doing some float math with integers...
        real(p2) :: ndist ! distance from wall node to cell center
        real(p2) :: xc, yc, zc, xn, yn, zn
        integer  :: n1,n2,n3,n4
        
        real(p2), dimension(3) :: fpoint
        
        continue

        ! Initialize the cell wall distance array
        allocate(cell_wall_distance(ncells))

        ! Define the pointer root_box
        allocate(root_box)
        allocate(gnode_to_wnode(nnodes))
        allocate(nf(nnodes))
        gnode_to_wnode = 0
        nf = 0  
        
        ! Build the array of interior nodes.
        is_wall = .false.
        bloop1 : do ib = 1,nb
            if (trim(bc_type(ib)) /= 'slip_wall' .and. trim(bc_type(ib)) /= 'no_slip_wall') cycle bloop1
            
            do inode = 1,bound(ib)%nbfaces
                do iface = 2,bound(ib)%bfaces(1,inode) + 1
                    ni = bound(ib)%bfaces(iface,inode)
                    is_wall(ni) = .true.
                    nf(ni) = nf(ni) + 1
                end do
            end do
        end do bloop1

        nwall_nodes = count(is_wall)

        ! Return early in the case where are no wall cells.
        if ( nwall_nodes == 0 ) then
            cell_wall_distance = my_huge
            return
        endif


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
        allocate(wn_to_wf(nwall_nodes))
        ! allocate(wn_sorty(nwall_nodes))
        ! allocate(wn_sortz(nwall_nodes))

        nwall_nodes = 0
        is_wall = .false.
        bloop2 : do ib = 1,nb
            if (trim(bc_type(ib)) /= 'slip_wall' .and. trim(bc_type(ib)) /= 'no_slip_wall') cycle bloop2
            
            do iface = 1,bound(ib)%nbfaces
                add_node_loop : do inode = 2,bound(ib)%bfaces(1,iface) + 1
                    ni = bound(ib)%bfaces(inode,iface)
                    if (.not.is_wall(ni)) then
                        is_wall(ni) = .true.
                        nwall_nodes = nwall_nodes + 1
                        wall_nodes(nwall_nodes) = ni
                        gnode_to_wnode(ni) = nwall_nodes
                        wnx(nwall_nodes)        = x(ni)
                        wny(nwall_nodes)        = y(ni)
                        wnz(nwall_nodes)        = z(ni)
                        
                        wn_to_wf(nwall_nodes)%nfaces = nf(ni)
                        allocate(wn_to_wf(nwall_nodes)%bface(nf(ni)))
                        allocate(wn_to_wf(nwall_nodes)%bound(nf(ni)))
                        nf(ni) = 0
                        wn_to_wf(nwall_nodes)%bound = 0
                        wn_to_wf(nwall_nodes)%bface = 0
                        wn_to_wf(nwall_nodes)%ni = ni
                    endif
                    nf(ni) = nf(ni) + 1
                    wni = gnode_to_wnode(ni)
                    wn_to_wf(wni)%bound(nf(ni)) = ib
                    wn_to_wf(wni)%bface(nf(ni)) = iface
                end do add_node_loop
            end do
        end do bloop2

        deallocate(nf)

        ! Set the split direction for the first bounding box
        call sort_longest(nwall_nodes,wnx,wny,wnz,wall_nodes,root_box,wns)

        call construct_bounding_box(nwall_nodes,wall_nodes,root_box)

        allocate(bbox_leafs(nleafs))

        nleafs = 0
        call extract_leafs(nleafs, root_box, bbox_leafs)

        deallocate(root_box)

        ! now we compute the distance to each box
        allocate(interior_cells(nleafs))
        allocate(icell_box_dist(nleafs))
        allocate(tmpinterior_cells(nleafs))
        allocate(tmpicell_box_dist(nleafs))

        cloop : do icell = 1,ncells
            do ibox = 1,nleafs
                tmpinterior_cells(ibox) = ibox
                tmpicell_box_dist(ibox) = distance_to_block(bbox_leafs(ibox)%xmin, bbox_leafs(ibox)%xmax, &
                                                               bbox_leafs(ibox)%ymin, bbox_leafs(ibox)%ymax, & 
                                                               bbox_leafs(ibox)%zmin, bbox_leafs(ibox)%zmax, &
                                                               cell(icell)%xc, cell(icell)%yc, cell(icell)%zc )
            end do
            ! fortran doesn't like it when you alias variables in subroutines.  I suppose I could rewrite the heap sort algorithm
            ! with an inout var but I don't want to right now...
            call heap_sort_index(tmpicell_box_dist,tmpinterior_cells,nleafs,interior_cells(:)) 
            do ibox = 1,nleafs
                ib = interior_cells(ibox)
                icell_box_dist(ibox) = tmpicell_box_dist(ib)
            end do

            ! Now that we have presorted the boxes by distance we compute the wall distance to the nodes in the closest box.  If the
            ! distance is within a threshold we compute the closest point on each surrounding face.  Since this is expensive we only
            ! do it for smaller wall distances.  Inverse wall distance is the only therm that shows up in any turbulence model which
            ! means past a certain point, "good enough" is fine.  Then we move to the next bounding box.  If the closest point is 
            ! greater than the wall distance we know we have the closest wall point.  If not we repeat
            cell_wall_distance(icell) = my_huge
            bloop3 : do ibox = 1,nleafs
                if (icell_box_dist(ibox) > cell_wall_distance(icell)) exit bloop3
                bxi = interior_cells(ibox)
                do inode = 1,bbox_leafs(bxi)%nwnodes
                    ni = bbox_leafs(bxi)%wnodes(inode)
                    xc = cell(icell)%xc
                    yc = cell(icell)%yc
                    zc = cell(icell)%zc
                    xn = x(ni)
                    yn = y(ni)
                    zn = z(ni)
                    ndist = sqrt( (xc-xn)**2 + (yc-yn)**2 + (zc-zn)**2)
                    if ( ndist < cell_wall_distance(icell) ) then
                        cell_wall_distance(icell) = ndist
                        closest_node = ni
                    endif
                end do
            end do bloop3

            if (cell_wall_distance(icell) < EXACT_DISTANCE_THRESHOLD) then
                ! compute a more exact wall distance
                wni = gnode_to_wnode(closest_node)
                do iface = 1,wn_to_wf(wni)%nfaces
                    fi = wn_to_wf(wni)%bface(iface)
                    ib = wn_to_wf(wni)%bound(iface)
                    select case(bound(ib)%bfaces(1,fi))
                    case(3) ! triangle
                        xc = cell(icell)%xc
                        yc = cell(icell)%yc
                        zc = cell(icell)%zc
                        n1 = bound(ib)%bfaces(2,fi)
                        n2 = bound(ib)%bfaces(3,fi)
                        n3 = bound(ib)%bfaces(4,fi)
                        fpoint = closestPointTriangle((/xc,yc,zc/),(/x(n1),y(n1),z(n1)/), &
                                                                   (/x(n2),y(n2),z(n2)/), &
                                                                   (/x(n3),y(n3),z(n3)/))
                    case(4) ! quad
                        n1 = bound(ib)%bfaces(2,fi)
                        n2 = bound(ib)%bfaces(3,fi)
                        n3 = bound(ib)%bfaces(4,fi)
                        n4 = bound(ib)%bfaces(5,fi)
                        fpoint = closestPointQuad((/xc,yc,zc/),(/x(n1),y(n1),z(n1)/), &
                                                               (/x(n2),y(n2),z(n2)/), &
                                                               (/x(n3),y(n3),z(n3)/), &
                                                               (/x(n4),y(n4),z(n4)/) )
                    case default 
                        write(*,*) "Error in the number of face sides"
                        write(*,*) "Stop. compute_wall_distance() wall_distance.f90"
                        stop
                    end select
                    ndist = sqrt( (xc-fpoint(1))**2 + (yc-fpoint(2))**2 + (zc-fpoint(3))**2)
                    cell_wall_distance(icell) = min(cell_wall_distance(icell),ndist)
                    
                end do
            end if

        end do cloop

        deallocate(wall_nodes)
        deallocate(wnx)
        deallocate(wny)
        deallocate(wnz)
        deallocate(wns)
        do inode = 1,nwall_nodes
            deallocate(wn_to_wf(inode)%bface)
            deallocate(wn_to_wf(inode)%bound)
        end do
        deallocate(wn_to_wf)


        deallocate(bbox_leafs)

        deallocate(interior_cells)
        deallocate(icell_box_dist)
        deallocate(tmpinterior_cells)
        deallocate(tmpicell_box_dist)
        deallocate(gnode_to_wnode)
    end subroutine compute_wall_distance

    pure function closestPointQuad(p,a,b,c,d) result(dist)

        implicit none
            
        real(p2), dimension(3), intent(in) :: p, a, b, c, d

        real(p2), dimension(3)             :: dist

        real(p2), dimension(3)  :: c1, c2
        real(p2)                :: d1, d2

        c1 = closestPointTriangle(p,a,b,c)
        d1 = sqrt( (p(1)-c1(1))**2 + (p(2)-c1(2))**2 + (p(3)-c1(3))**2 )

        c2 = closestPointTriangle(p,a,c,d)
        d2 = sqrt( (p(1)-c2(1))**2 + (p(2)-c2(2))**2 + (p(3)-c2(3))**2 )
        
        if (d1 < d2) then 
            dist = c1
        else
            dist = c2
        endif

    end function closestPointQuad

    pure function closestPointTriangle(p,a,b,c) result(dist)
        ! This function is "borrowed" (stolen) from Intel's open source embree graphics library.
        ! why do the work again when someone already did it for me?
        ! what's the saying? Good musicians borrow, great musicians steal...
        ! https://github.com/RenderKit/embree/blob/master/tutorials/common/math/closest_point.h

        implicit none
        
        real(p2), dimension(3), intent(in) :: p, a, b, c

        real(p2), dimension(3)             :: dist

        real(p2), dimension(3) :: ab, ac, ap, bp, cp
        real(p2)               :: d1,d2,d3,d4,d5,d6
        real(p2)               :: vc, v, vb, va, w
        real(p2)               :: denom

        ab = b - a
        ac = c - a
        ap = p - a

        d1 = dot_product(ab,ap)
        d2 = dot_product(ac,ap)
        if (d1 <= 0.0_p2 .and. d2 <= 0.0_p2) then
            dist = a
            return
        end if

        bp = p - b
        d3 = dot_product(ab,bp)
        d4 = dot_product(ac,bp)
        if (d3 >= 0.0_p2 .and. d4 <= d3) then
            dist = b
            return
        end if

        cp = p - c
        d5 = dot_product(ab,cp)
        d6 = dot_product(ac,cp)
        if (d6 >= 0.0_p2 .and. d5 <= d6) then
            dist = c
            return
        end if

        vc = d1 * d4 - d3 * d2
        if (vc <= 0.0_p2 .and. d1 >= 0.0_p2 .and. d3 <= 0.0_p2) then
            v = d1 / (d1 - d3)
            dist = a + v * ab
            return
        end if

        vb = d5 * d2 - d1 * d6
        if (vb <= 0.0_p2 .and. d2 >= 0.0_p2 .and. d6 <= 0.0_p2) then
            V = d2 / (d2 - d6)
            dist = a + v * ac
            return
        end if

        va = d3 * d6 - d5 * d4
        if (va <= 0.0_p2 .and. (d4 - d3) >= 0.0_p2 .and. (d5 - d6) >= 0.0_p2) then
            v = (d4 - d3) / ((d4 - d3) + (d5 - d6))
            dist = b + v * (c - b)
        endif

        denom = 1.0_p2 / (va + vb + vc)
        v = vb * denom
        w = vc * denom
        dist = a + v * ab + w * ac
        
    end function closestPointTriangle

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

        integer :: split, nnms
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
            if (associated(leaf_array(nlf)%branch1)) deallocate(leaf_array(nlf)%branch1)
            if (associated(leaf_array(nlf)%branch2)) deallocate(leaf_array(nlf)%branch2)
        else
            call extract_leafs(nlf, root_box%branch1, leaf_array)
            call extract_leafs(nlf, root_box%branch2, leaf_array)
        end if

        if (associated(root_box%root)) nullify(root_box%root)
        if (associated(root_box%branch1)) deallocate(root_box%branch1)
        if (associated(root_box%branch2)) deallocate(root_box%branch2)
        if (allocated( root_box%wnodes )) deallocate(root_box%wnodes)

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