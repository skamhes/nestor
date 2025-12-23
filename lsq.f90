module least_squares
    
    ! Module for computing the gradient using least squares.  Available options
    ! 
    ! Weighted vertex gradient: https://doi.org/10.48550/arXiv.1702.04518
    ! Computes the gradient at each vertex and then uses them to find an average gradient at the cell.  Since the field value
    ! is not know at the vertex it can also be solved for, resulting in 4 unknowns for each node.  Currently that value is not used.
    ! At boundaries the treatment is a little different.  To ensure a defined (or over-defined since it is a LSQ method) system
    ! ghost cells are used to provide additional rows to the linear system.  Additionally, depending on the BC type, some field
    ! values may be known a priori (such as zero velocity at the walls).  In this instence we move the vertex value to the RHS and 
    ! use a different set of LSQ coefficients to compute the gradient.  This will result in a slightly larger memory footprint but 
    ! presumabley should result in increased robustness.

    use common , only : p2

    implicit none

    public

    type lsq_vertex_type
        integer                           ::   ncells_lsq  ! number of cells attached to lsq vertex
        integer , dimension(:)  , pointer ::     cell_lsq  ! lsq cell
        integer , dimension(:)  , pointer ::       ib_lsq  ! boundary number for each attached cell (0 = internal cell)
        real(p2), dimension(:)  , pointer ::          cq4  ! LSQ coefficient for q at vertex (4 unknowns)
        real(p2), dimension(:)  , pointer ::          cx4  ! LSQ coefficient for x-derivative (4 unknowns)
        real(p2), dimension(:)  , pointer ::          cy4  ! LSQ coefficient for y-derivative (4 unknowns) 
        real(p2), dimension(:)  , pointer ::          cz4  ! LSQ coefficient for z-derivative (4 unknowns) 
        real(p2), dimension(:)  , pointer ::          cx3  ! LSQ coefficient for x-derivative (4 unknowns) 
        real(p2), dimension(:)  , pointer ::          cy3  ! LSQ coefficient for y-derivative (4 unknowns) 
        real(p2), dimension(:)  , pointer ::          cz3  ! LSQ coefficient for z-derivative (4 unknowns) 
        integer                           ::        btype  ! See below for values. 0 = internal cell
    end type lsq_vertex_type

    !Cell data array in the custom data type.
    type(lsq_vertex_type), dimension(:), pointer :: lsqv  !cell-centered LSQ array

    type lsq_cell_type
        integer                                 ::    n_nnghbrs  ! number of cells attached to lsq vertex
        integer,  dimension(:)  , pointer       ::    nghbr_lsq  ! list of neighbor cells  
        real(p2), dimension(:,:,:), pointer     ::           cf  ! LSQ coefficient for x,y,&z-derivative (3 unknowns) 
        real(p2), dimension(:,:), pointer     ::           cx  ! LSQ coefficient for x-derivative (3 unknowns) 
        real(p2), dimension(:,:), pointer     ::           cy  ! LSQ coefficient for y-derivative (3 unknowns) 
        real(p2), dimension(:,:), pointer     ::           cz  ! LSQ coefficient for z-derivative (3 unknowns) 
        integer                                 ::          nbf  ! number of boundary faces attached to the cell
        integer,  dimension(:,:), pointer       ::       gcells  ! list of ghost cells (ibcell,ib) length 2xnbf
        real(p2), dimension(:,:,:), pointer     ::          gcf  ! LSQ coefficient for x,y,&z-derivative (3 unknowns) 
        real(p2), dimension(:,:), pointer     ::          gcx  ! LSQ coefficient for x-derivative (3 unknowns) 
        real(p2), dimension(:,:), pointer     ::          gcy  ! LSQ coefficient for y-derivative (3 unknowns) 
        real(p2), dimension(:,:), pointer     ::          gcz  ! LSQ coefficient for z-derivative (3 unknowns) 
    end type lsq_cell_type

    !Cell data array in the custom data type.
    type(lsq_cell_type), dimension(:), pointer :: lsqc  !cell-centered LSQ array

    type glsq_type 
        type(lsq_cell_type), dimension(:), pointer :: lsq ! this type is just variable length arrays of lsq cell types
    end type glsq_type

    type(glsq_type), dimension(:), pointer :: lsqg


    ! Boundary int rankings
    integer, parameter :: NO_SLIP_WALL = 99
    integer, parameter :: SLIP_WALL = 98
    integer, parameter :: FREE_STREAM = 51
    ! IF BTYPE > 50 THEN NODE_VALUE = KNOWN
    integer, parameter :: PRESSURE_OUTLET = 20
    integer, parameter :: INTERNAL = 0

    public construct_lsq_stencil

    contains

    subroutine construct_lsq_stencil

        ! use common , only : p2

        use utils , only : ilsq_stencil, LSQ_STENCIL_WVERTEX, LSQ_STENCIL_NN

        implicit none

        select case(ilsq_stencil)
        case(LSQ_STENCIL_WVERTEX)
            call construct_wvertex_stencil
            call compute_vertex_coefficients
        case(LSQ_STENCIL_NN)
            call construct_nn_stencil
            call compute_cell_coefficients
        case default
            write(*,*) "Unsupported LSQ Stencil"
            stop
        end select



    end subroutine construct_lsq_stencil

    subroutine construct_wvertex_stencil

        use grid , only : nnodes, ncells, bound, nb, cell, bc_type

        implicit none

        !To store the list of cells around each node.
        type node_type
            integer                        :: nic ! number of internal cells
            integer, dimension(:), pointer :: c   ! internal cells
            integer                        :: nbc ! number of boundary cells (0 = internal)
            integer, dimension(:), pointer :: bc  ! boundary cells
        end type node_type

        type(node_type), dimension(:), pointer :: node

        integer :: i, j, k, ib
        integer :: vk, fj !, ck

        write(*,*)
        write(*,*) " --------------------------------------------------"
        write(*,*) " Constructing vertex neighbors... "
        write(*,*)

        allocate( node(nnodes) )
        allocate( lsqv(nnodes ) )

        do i = 1,nnodes
            node(i)%nic = 0
            node(i)%nbc = 0
        end do

        ! Count the number of internal cells
        do i = 1, ncells
            do k = 1, cell(i)%nvtx
                vk           = cell(i)%vtx(k)
                node(vk)%nic = node(vk)%nic + 1
            end do
        end do

        ! Count the number of boundary nodes
        do ib = 1,nb
            do j = 1,bound(ib)%nbfaces
                fj           = bound(ib)%bfaces(1,j) ! (/#_of_nodes v1:v_end/)
                do k = 2,(fj+1) 
                    vk           = bound(ib)%bfaces(k,j)
                    node(vk)%nbc = node(vk)%nbc + 1
                end do
            end do
        end do

        ! Allocate arrays in the lsq struct array
        do i = 1,nnodes
            lsqv(i)%ncells_lsq = node(i)%nic + node(i)%nbc
            allocate( lsqv(i)%cell_lsq( lsqv(i)%ncells_lsq  ) )
            allocate( lsqv(i)%ib_lsq( lsqv(i)%ncells_lsq ) )
            lsqv(i)%btype = 0 ! stays zero for internal cells

            node(i)%nic = 0
            node(i)%nbc = 0
        end do

        ! Add internal cells to the array
        do i = 1,ncells
            do k = 1, cell(i)%nvtx
                vk           = cell(i)%vtx(k)
                node(vk)%nic = node(vk)%nic + 1
                lsqv(vk)%cell_lsq(node(vk)%nic) = i
                lsqv(vk)%ib_lsq(node(vk)%nic) = INTERNAL ! Initialize the cell as internal
            end do
        end do

        ! Add boundary cells to the array
        do ib = 1,nb
            do j = 1,bound(ib)%nbfaces
                fj           = bound(ib)%bfaces(1,j) ! (/#_of_nodes v1:v_end/)
                do k = 2,(fj+1) 
                    vk           = bound(ib)%bfaces(k,j)
                    node(vk)%nbc = node(vk)%nbc + 1
                    lsqv(vk)%cell_lsq(node(vk)%nic + node(vk)%nbc) = j
                    ! We will need face data in some cases so we will use bface# in place of cell
                    ! In the future I can sort these so the boundary faces match with their attached cell
                    ! In theory this would allow me to skip some assignmets.  Will do another time.
                    
                    ! Give boundary of attached bcell
                    lsqv(vk)%ib_lsq(node(vk)%nic + node(vk)%nbc) = ib ! internal cell

                    ! Assign the boundary condition based on above hierarchy
                    select case(trim(bc_type(ib))) 
                        case('freestream')
                            lsqv(vk)%btype = max(lsqv(vk)%btype,FREE_STREAM)
                        case('symmetry')
                            ! symmetry is treated numerically like a slip_wall
                            lsqv(vk)%btype = max(lsqv(vk)%btype,SLIP_WALL)
                        case('slip_wall')
                            lsqv(vk)%btype = max(lsqv(vk)%btype,SLIP_WALL)
                        case('no_slip_wall')
                            lsqv(vk)%btype = max(lsqv(vk)%btype,NO_SLIP_WALL)
                        case('outflow_subsonic')
                            lsqv(vk)%btype = max(lsqv(vk)%btype,PRESSURE_OUTLET)
                        case default
                            write(*,*) "Boundary condition=",trim(bc_type(ib)),"  not implemented."
                            stop
                    end select
                end do
            end do
        end do

        ! We no longer need the node array
        deallocate(node)

        ! Allocate the cell weight values
        do i = 1,nnodes
            allocate( lsqv(i)%cx4(lsqv(i)%ncells_lsq) )
            allocate( lsqv(i)%cy4(lsqv(i)%ncells_lsq) )
            allocate( lsqv(i)%cz4(lsqv(i)%ncells_lsq) )
            allocate( lsqv(i)%cq4(lsqv(i)%ncells_lsq) )
            if ( lsqv(i)%btype > 0 ) then
                allocate( lsqv(i)%cx3(lsqv(i)%ncells_lsq) )
                allocate( lsqv(i)%cy3(lsqv(i)%ncells_lsq) )
                allocate( lsqv(i)%cz3(lsqv(i)%ncells_lsq) )
            endif
        end do
        ! I want the option to prescribe the node value at boundaries where it's known.
        ! It seems like this should provide an increased robustness(?).  Intuition says I can't
        ! just use the calculated weights for cy, cx, and cz, but I'll have to investigate...
        
    end subroutine construct_wvertex_stencil

    subroutine construct_nn_stencil

        use common , only : p2

        use grid , only : cell, x, y, z, nnodes, bound, ncells, node_type, nb

        use sort_routines , only : queued_natural_merge_sort

        use utils , only : iturb_type, TURB_INVISCID

        use solution , only : nlsq

        implicit none

        type(node_type), dimension(nnodes) :: node

        type bnode_type
            type(node_type), dimension(:), pointer :: node
        end type bnode_type

        type(bnode_type), dimension(nb) :: bnode
    
        type nn_type
            integer                             :: n_nnghbr
            integer, dimension(:)  , pointer    :: nnghbr
        end type nn_type

        type(nn_type), dimension(ncells) :: c2nn ! Array that lists the node neighbors of each cell
        integer, dimension(ncells) :: c2bf ! pointer from cell to bface

        integer, dimension(:), pointer :: scratch_nghbrs
        integer                        :: dupn_nnghbr

        integer, dimension(:), pointer :: runpointer
        
        integer :: inode, icell, ib
        integer :: jcell, cj
        integer :: ni, cnvtx, ci, cn, nbfn

        integer :: start, end

        if (iturb_type == TURB_INVISCID) then
            nlsq = 1
        else
            nlsq = 2 
        endif

        allocate(lsqc(ncells))
        c2bf = 0

        do inode = 1,nnodes
            node(inode)%nc = 0
        end do

        do icell = 1,ncells
            do inode = 1,cell(icell)%nvtx
                ni = cell(icell)%vtx(inode)
                node(ni)%nc = node(ni)%nc + 1
            end do
        end do

        do inode = 1,nnodes
            allocate(node(inode)%c(node(inode)%nc))
            node(inode)%nc = 0
        end do

        do icell = 1,ncells
            do inode = 1,cell(icell)%nvtx
                ni = cell(icell)%vtx(inode)
                node(ni)%nc = node(ni)%nc + 1
                node(ni)%c(node(ni)%nc)  = icell
            end do
        end do

        ! At this point we have a vector of all the nodes with a (sorted) list and # of attached cells
        
        allocate(scratch_nghbrs(8))
        allocate(runpointer(8))
        
        do icell = 1,ncells

            cnvtx = cell(icell)%nvtx

            ! Count the number of node neighbors including duplicates
            dupn_nnghbr = 0            
            do inode = 1,cnvtx
                ni = cell(icell)%vtx(inode)
                dupn_nnghbr = dupn_nnghbr + node(ni)%nc 
            end do

            ! it's ok if the scratch vector is too long but we don't want an overflow
            if (size(scratch_nghbrs) < dupn_nnghbr) then
                deallocate(scratch_nghbrs)
                allocate(scratch_nghbrs(dupn_nnghbr))
            endif

            ! Add the node neighbors to the scratch vetor
            dupn_nnghbr = 0
            if (cnvtx >= size(runpointer)) then ! this should be cnvtx + 1 > rp but for ints n+1>x <=> n>=x
                deallocate(runpointer)
                allocate(runpointer(cnvtx+1))
            endif
            runpointer(1) = 1
            do inode = 1,cnvtx
                ni = cell(icell)%vtx(inode)
                start = dupn_nnghbr + 1
                end   = dupn_nnghbr + node(ni)%nc
                scratch_nghbrs(start:end) = node(ni)%c(:)
                dupn_nnghbr = dupn_nnghbr + node(ni)%nc 
                runpointer(inode + 1) = end + 1 ! this will be used by the reduction algorithm
            end do

            call queued_natural_merge_sort(dupn_nnghbr,cnvtx,runpointer, &
                                            scratch_nghbrs, icell,c2nn(icell)%n_nnghbr, c2nn(icell)%nnghbr)
            
            lsqc(icell)%n_nnghbrs = c2nn(icell)%n_nnghbr ! Not sure why I included this step...
            allocate(lsqc(icell)%nghbr_lsq(c2nn(icell)%n_nnghbr))
            lsqc(icell)%nghbr_lsq = c2nn(icell)%nnghbr
            
            allocate(lsqc(icell)%cx(c2nn(icell)%n_nnghbr,nlsq))
            allocate(lsqc(icell)%cy(c2nn(icell)%n_nnghbr,nlsq))
            allocate(lsqc(icell)%cz(c2nn(icell)%n_nnghbr,nlsq))
            allocate(lsqc(icell)%cf(3,c2nn(icell)%n_nnghbr,nlsq))

            lsqc(icell)%nbf = 0
        end do

        

        ! Now we need to add the ghost cells
        allocate(lsqg(nb))
        do ib = 1,nb
            allocate(lsqg(ib)%lsq(bound(ib)%nbfaces))
            ! Create a pointer that maps a cell to its bound(ib)%bcell() index.
            do icell = 1,bound(ib)%nbfaces
                ci = bound(ib)%bcell(icell)
                c2bf(ci) = icell ! pointer maps cell index to bcell index
            end do
            do icell = 1,bound(ib)%nbfaces ! loop through the bcells
                ci = bound(ib)%bcell(icell)
                if (size(scratch_nghbrs) < lsqc(ci)%n_nnghbrs) then
                    deallocate(scratch_nghbrs)
                    allocate(scratch_nghbrs(lsqc(ci)%n_nnghbrs))
                endif
                lsqg(ib)%lsq(icell)%n_nnghbrs = 1 ! the ghost attached to bface(ci) counts as 1
                scratch_nghbrs(1) = icell
                do jcell = 1,lsqc(ci)%n_nnghbrs ! loop through the internal neighbors to ci which we already know
                    cj = lsqc(ci)%nghbr_lsq(jcell)
                    if (c2bf(cj) > 0) then ! if the map is nonzero for that cell it is also on the given boundary
                        lsqg(ib)%lsq(icell)%n_nnghbrs = lsqg(ib)%lsq(icell)%n_nnghbrs + 1 ! increment
                        scratch_nghbrs(lsqg(ib)%lsq(icell)%n_nnghbrs) = c2bf(cj) ! add it to the list of node neighbors
                        ! because we're looping over the merged list we don't need to merge again
                    endif
                end do
                allocate(lsqg(ib)%lsq(icell)%nghbr_lsq(lsqg(ib)%lsq(icell)%n_nnghbrs))
                lsqg(ib)%lsq(icell)%nghbr_lsq = scratch_nghbrs(1:lsqg(ib)%lsq(icell)%n_nnghbrs)
                lsqc(ci)%nbf = lsqc(ci)%nbf + lsqg(ib)%lsq(icell)%n_nnghbrs
            end do
            do icell = 1,bound(ib)%nbfaces
                ! reset the c2bf array to 0
                ci = bound(ib)%bcell(icell)
                c2bf(ci) = 0
            end do
        end do

        do ib = 1,nb
            do icell = 1,bound(ib)%nbfaces
                ci = bound(ib)%bcell(icell)
                if (.not.associated(lsqc(ci)%gcells)) then 
                    allocate(lsqc(ci)%gcells(2,lsqc(ci)%nbf))
                    allocate(lsqc(ci)%gcx(lsqc(ci)%nbf,nlsq))
                    allocate(lsqc(ci)%gcy(lsqc(ci)%nbf,nlsq))
                    allocate(lsqc(ci)%gcz(lsqc(ci)%nbf,nlsq))
                    allocate(lsqc(ci)%gcf(3,lsqc(ci)%nbf,nlsq))
                    lsqc(ci)%nbf = 0
                end if
                
                do jcell = 1,lsqg(ib)%lsq(icell)%n_nnghbrs
                    lsqc(ci)%nbf = lsqc(ci)%nbf + 1
                    lsqc(ci)%gcells(:,lsqc(ci)%nbf) = (/ lsqg(ib)%lsq(icell)%nghbr_lsq(jcell) , ib /)    
                end do
                deallocate(lsqg(ib)%lsq(icell)%nghbr_lsq)
            end do
            deallocate(lsqg(ib)%lsq)
        end do
        deallocate(lsqg)
        deallocate(scratch_nghbrs)
    end subroutine construct_nn_stencil

    subroutine compute_vertex_coefficients

        use grid , only : cell, x, y, z, nnodes, bound, ncells

        use common , only : p2, zero, one, two

        use direct_solve , only : qr_factorization

        implicit none

        real(p2) :: maxdx, maxdy, maxdz
        real(p2) :: lsq_weight_invdis_power
        integer                           :: m!, n             !Size of LSQ matrix: A(m,n).
        real(p2), pointer, dimension(:,:) :: a3, a4           !LSQ matrix: A(m,n).
        real(p2), pointer, dimension(:,:) :: rinvqt3, rinvqt4 !Pseudo inverse R^{-1}*Q^T
        integer                           :: connect_cell, connect_bface
        real(p2), pointer, dimension(:,:) :: cell_grad
        
        integer :: i, k, ib
        ! integer :: ck, fk
        
        real(p2) :: dx, dy, dz
        real(p2) :: cgx, cgy, cgz ! ghost cell "center"
        real(p2) :: weight_k
        logical  :: verification_error
        real(p2) :: wx, wy, wz, wq
        real(p2) :: xi, yi, zi
        real(p2) :: xk, yk, zk

        real(p2), dimension(3) :: maxDeltasNZ

        ! Debug:
        integer :: unknowns_ = 3

        write(*,*)
        write(*,*) "--------------------------------------------------"
        write(*,*) " Computing LSQ coefficients... "
        write(*,*)

        maxdx = zero
        maxdy = zero
        maxdz = zero

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! The power to the inverse distance weight. The value 0.0 is used to avoid
        ! instability known for Euler solvers. So, this is the unweighted LSQ gradient.
        ! More accurate gradients are obtained with 1.0, and such can be used for the
        ! viscous terms and source terms in turbulence models.
        lsq_weight_invdis_power = 0

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! Compute the LSQ coefficients (cq,cx,cy,cz) at all nodes.
        node_loop : do i=1,nnodes
            m = lsqv(i)%ncells_lsq ! # of connected cells
            allocate( a3(m,3) )
            allocate( a4(m,4) )
            allocate( rinvqt3(3,m) )
            allocate( rinvqt4(4,m) )

            connect_loop : do k = 1,m
                if ( lsqv(i)%ib_lsq(k) == INTERNAL ) then ! Internal ib = 0
                    connect_cell = lsqv(i)%cell_lsq(k)
                    dx = cell(connect_cell)%xc - x(i)
                    dy = cell(connect_cell)%yc - y(i)
                    dz = cell(connect_cell)%zc - z(i)
                else
                    connect_bface = lsqv(i)%cell_lsq(k)
                    ib            = lsqv(i)%ib_lsq(k)
                    connect_cell  = bound(ib)%bcell(connect_bface)
                    ! briefly use the dx, dy, znd dz for face_center - cell_center
                    dx = bound(ib)%bface_center(1,connect_bface) - cell(connect_cell)%xc
                    dy = bound(ib)%bface_center(2,connect_bface) - cell(connect_cell)%yc
                    dz = bound(ib)%bface_center(3,connect_bface) - cell(connect_cell)%zc
                    cgx = bound(ib)%bface_center(1,connect_bface) + dx
                    cgy = bound(ib)%bface_center(2,connect_bface) + dy
                    cgz = bound(ib)%bface_center(3,connect_bface) + dz
                    dx = cgx - x(i)
                    dy = cgy - y(i)
                    dz = cgz - z(i)
                endif
                weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                ! 3 unknowns
                a3(k,1) = weight_k * dx
                a3(k,2) = weight_k * dy
                a3(k,3) = weight_k * dz
                ! 4 unknowns
                a4(k,1) = weight_k !* one
                a4(k,2) = weight_k * dx
                a4(k,3) = weight_k * dy
                a4(k,4) = weight_k * dz
                ! 2D check
                maxdx  = max(abs(dx),maxdx)
                maxdy  = max(abs(dy),maxdy)
                maxdz  = max(abs(dz),maxdz)
            end do connect_loop
            !-------------------------------------------------------
            ! Perform QR factorization and compute R^{-1}*Q^T from A(m,n).
            call qr_factorization(a4,rinvqt4,m,4)
            call qr_factorization(a3,rinvqt3,m,3)

            !-------------------------------------------------------
            ! Compute and store the LSQ coefficients: R^{-1}*Q^T*w.
            !
            ! (wx,wy,wz) = R^{-1}*Q^T*RHS
            !            = sum_k (cx,cy,cz)*(wk-wi).
            connect_loop2 : do k = 1,m
                if ( lsqv(i)%ib_lsq(k) == INTERNAL ) then ! Internal ib = 0
                    connect_cell = lsqv(i)%cell_lsq(k)
                    dx = cell(connect_cell)%xc - x(i)
                    dy = cell(connect_cell)%yc - y(i)
                    dz = cell(connect_cell)%zc - z(i)
                else
                    connect_bface = lsqv(i)%cell_lsq(k)
                    ib            = lsqv(i)%ib_lsq(k)
                    connect_cell  = bound(ib)%bcell(connect_bface)
                    ! briefly use the dx, dy, znd dz for face_center - cell_center
                    dx = bound(ib)%bface_center(1,connect_bface) - cell(connect_cell)%xc
                    dy = bound(ib)%bface_center(2,connect_bface) - cell(connect_cell)%yc
                    dz = bound(ib)%bface_center(3,connect_bface) - cell(connect_cell)%zc
                    cgx = bound(ib)%bface_center(1,connect_bface) + dx
                    cgy = bound(ib)%bface_center(2,connect_bface) + dy
                    cgz = bound(ib)%bface_center(3,connect_bface) + dz
                    dx = cgx - x(i)
                    dy = cgy - y(i)
                    dz = cgz - z(i)
                endif
                weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                ! 3 unknowns
                lsqv(i)%cx3(k) = rinvqt3(1,k) * weight_k
                lsqv(i)%cy3(k) = rinvqt3(2,k) * weight_k
                lsqv(i)%cz3(k) = rinvqt3(3,k) * weight_k
                ! 4 unknowns
                lsqv(i)%cq4(k) = rinvqt4(1,k) * weight_k
                lsqv(i)%cx4(k) = rinvqt4(2,k) * weight_k
                lsqv(i)%cy4(k) = rinvqt4(3,k) * weight_k
                lsqv(i)%cz4(k) = rinvqt4(4,k) * weight_k
            end do connect_loop2

            deallocate(a3,rinvqt3)
            deallocate(a4,rinvqt4)

        end do node_loop

        ! Verification
        ! Compute the gradient of w = 2*x+y+4*z to se if we get wx = 2, wy = 1, and wz = 4 correctly
        verification_error = .false.

        ! Initialize cell centered gradient
        allocate(cell_grad(3,ncells))
        cell_grad = zero

        ! First calculate the gradient at each node.
        do i = 1, nnodes
            wx = zero
            wy = zero
            wz = zero
            wq = zero
            ! (xi,yi,zi) to be used to compute the function 2*x+y+4z at i
            xi = x(i)
            yi = y(i)
            zi = z(i)

            ! Loop over neighbor cells
            do k = 1,lsqv(i)%ncells_lsq
                if ( lsqv(i)%ib_lsq(k) == INTERNAL ) then
                    connect_cell = lsqv(i)%cell_lsq(k)
                    xk = cell(connect_cell)%xc
                    yk = cell(connect_cell)%yc
                    zk = cell(connect_cell)%zc
                else
                    connect_bface = lsqv(i)%cell_lsq(k) 
                    ib            = lsqv(i)%ib_lsq(k)
                    connect_cell  = bound(ib)%bcell(connect_bface)
                    ! We don't store these anywhere because we only need the ghost cell center for this check.
                    ! Ordinarilly we don't need it.
                    dx = bound(ib)%bface_center(1,connect_bface) - cell(connect_cell)%xc
                    dy = bound(ib)%bface_center(2,connect_bface) - cell(connect_cell)%yc
                    dz = bound(ib)%bface_center(3,connect_bface) - cell(connect_cell)%zc
                    xk = bound(ib)%bface_center(1,connect_bface) + dx ! cgx
                    yk = bound(ib)%bface_center(2,connect_bface) + dy ! cgy
                    zk = bound(ib)%bface_center(3,connect_bface) + dz ! cgz
                endif
                ! This is how we use the LSQ coefficients: accumulate cx*(wk-wi)
                ! and cy*(wk-wi) and cz*(wk-wi)
                if ( unknowns_ == 4 ) then
                    wx = wx + lsqv(i)%cx4(k)*( (2.0*xk+yk+4.0*zk) )
                    wy = wy + lsqv(i)%cy4(k)*( (2.0*xk+yk+4.0*zk) )
                    wz = wz + lsqv(i)%cz4(k)*( (2.0*xk+yk+4.0*zk) )
                    wq = wq + lsqv(i)%cq4(k)*( (2.0*xk+yk+4.0*zk) )
                    ! we don't need q
                else ! unknowns_ == 3
                    wx = wx + lsqv(i)%cx3(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yk+4.0*zi) )
                    wy = wy + lsqv(i)%cy3(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wz = wz + lsqv(i)%cz3(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                endif
            end do
            ! Loop through attached cells...again
            cell_grad_loop : do k = 1,lsqv(i)%ncells_lsq
                if ( lsqv(i)%ib_lsq(k) == INTERNAL ) then
                    connect_cell = lsqv(i)%cell_lsq(k)
                else
                    cycle cell_grad_loop
                    ! We should be able to cycle the node loop since boundary cells are added last.
                    ! But if I change the sorting for that array that would no longer be the case.  I'll leave it like this for now
                    ! as I don't mind the inefficiency for now, and I don't want to have to go bug hunting later...
                endif                
                cell_grad(1,connect_cell) = cell_grad(1,connect_cell) + wx
                cell_grad(2,connect_cell) = cell_grad(2,connect_cell) + wy
                cell_grad(3,connect_cell) = cell_grad(3,connect_cell) + wz
            end do cell_grad_loop
            maxDeltasNZ = zero
            if (maxdx > 0.001_p2) maxDeltasNZ(1) = one
            if (maxdy > 0.001_p2) maxDeltasNZ(2) = one
            if (maxdz > 0.001_p2) maxDeltasNZ(3) = one
            if ( maxDeltasNZ(1)*abs(wx-two) > 1.0e-06_p2 .or. &
                 maxDeltasNZ(2)*abs(wy-one) > 1.0e-06_p2 .or. &
                 maxDeltasNZ(3)*abs(wz-4.0_p2) > 1.0e-06_p2) then
                    write(*,*) " NODE: ", i
                    write(*,*) " wx = ", wx, " exact ux = 2.0"!,maxDeltasNZ(1)*abs(wx-two)
                    write(*,*) " wy = ", wy, " exact uy = 1.0"!,maxDeltasNZ(2)*abs(wy-one)
                    write(*,*) " wz = ", wz, " exact uz = 4.0"!, maxDeltasNZ(3)*abs(wz-4.0_p2),maxDeltasNZ(3)
                    if ( unknowns_ == 4 ) then
                        write(*,*) " wq = ", wq , " exact wq = ", (2.0*xi+yi+4.0*zi)
                    end if
                    verification_error = .true.
            end if
        end do
        if (verification_error) then

            write(*,*) " LSQ coefficients are not correct. See above. Stop."
            stop
         
        else
         
            write(*,*) " Verified: LSQ coefficients at node are exact for a linear function."
            write(*,*)
         
        endif

        do i = 1,ncells
            !                                 INT2REAL
            cell_grad(:,i) = cell_grad(:,i) / real(cell(i)%nvtx,p2)
            wx = cell_grad(1,i)
            wy = cell_grad(2,i)
            wz = cell_grad(3,i)
            maxDeltasNZ = zero
            if (maxdx > 0.001_p2) maxDeltasNZ(1) = one
            if (maxdy > 0.001_p2) maxDeltasNZ(2) = one
            if (maxdz > 0.001_p2) maxDeltasNZ(3) = one
            if ( maxDeltasNZ(1)*abs(wx-two) > 1.0e-06_p2 .or. &
                 maxDeltasNZ(2)*abs(wy-one) > 1.0e-06_p2 .or. &
                 maxDeltasNZ(3)*abs(wz-4.0_p2) > 1.0e-06_p2) then
                    write(*,*) " CELL: ", i
                    write(*,*) " wx = ", wx, " exact ux = 2.0"!,maxDeltasNZ(1)*abs(wx-two)
                    write(*,*) " wy = ", wy, " exact uy = 1.0"!,maxDeltasNZ(2)*abs(wy-one)
                    write(*,*) " wz = ", wz, " exact uz = 4.0"!, maxDeltasNZ(3)*abs(wz-4.0_p2),maxDeltasNZ(3)
                    verification_error = .true.
            end if
        end do

        if (verification_error) then

            write(*,*) " LSQ coefficients result in inaccurate cell values. See above. Stop."
            stop
         
        else
         
            write(*,*) " Verified: LSQ coefficients at cells are exact for a linear function."
         
        endif

        deallocate(cell_grad)
        
        write(*,*)
        write(*,*) " End of Computing LSQ coefficients... "
        write(*,*) "--------------------------------------------------"
        write(*,*)
        
        ! End of Compute the LSQ coefficients in all cells.
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
         
    end subroutine compute_vertex_coefficients

    subroutine compute_cell_coefficients
        
        use grid , only : cell, ncells, gcell

        use common , only : p2, zero, one, two, ix, iy, iz

        use direct_solve , only : qr_factorization

        use solution , only : ndim, nlsq

        implicit none

        real(p2) :: maxdx, maxdy, maxdz
        real(p2) :: lsq_weight_invdis_power
        integer                             :: m, n             !Size of LSQ matrix: A(m,n).
        real(p2), pointer, dimension(:,:,:) :: a                !LSQ matrix: A(m,n).
        real(p2), pointer, dimension(:,:,:) :: rinvqt           !Pseudo inverse R^{-1}*Q^T
        integer                             :: nghbr_cell
        
        integer :: k, ib, mweight
        integer :: icell
        integer :: ci
        
        real(p2) :: dx, dy, dz
        real(p2) :: weight_k
        logical  :: verification_error
        real(p2) :: wx, wy, wz
        real(p2) :: xi, yi, zi
        real(p2) :: xk, yk, zk

        real(p2), dimension(3) :: maxDeltasNZ

        write(*,*)
        write(*,*) "--------------------------------------------------"
        write(*,*) " Computing LSQ coefficients... "
        write(*,*)

        maxdx = zero
        maxdy = zero
        maxdz = zero

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! The power to the inverse distance weight. The value 0.0 is used to avoid
        ! instability known for Euler solvers. So, this is the unweighted LSQ gradient.
        ! More accurate gradients are obtained with 1.0, and such can be used for the
        ! viscous terms and source terms in turbulence models.
        ! lsq_weight_invdis_power = 0

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! Compute the LSQ coefficients (cq,cx,cy,cz) at all cells.
        cloop : do icell = 1,ncells
            m = lsqc(icell)%n_nnghbrs   ! # of neighbors
            n = ndim                    ! # of dimensions

            ! Allocate LSQ matrix and the pseudo inverse, R^{-1}*Q^T
            allocate(a(m + lsqc(icell)%nbf,n,nlsq)) ! note: it may produce some additional speed to switch the rows and columns here
            ! however a is a very small matrix and this subroutine is called once so we'll leave that for another day
            allocate(rinvqt(n,m + lsqc(icell)%nbf,nlsq))
            ! Initialize a
            a = zero
           
            !-------------------------------------------------------
            ! Build the weighted-LSQ matrix A(m,n).
            !
            !     weight_1 * [ (x1-xi)*wxi + (y1-yi)*wyi + (z1-zi)*wzi ] = weight_1 * [ w1 - wi ]
            !     weight_2 * [ (x2-xi)*wxi + (y2-yi)*wyi + (z2-zi)*wzi ] = weight_2 * [ w2 - wi ]
            !                 .
            !                 .
            !     weight_m * [ (xm-xi)*wxi + (ym-yi)*wyi + (zm-zi)*wzi ] = weight_2 * [ wm - wi ]
            nghbr_loop : do k = 1, lsqc(icell)%n_nnghbrs
                nghbr_cell = lsqc(icell)%nghbr_lsq(k) !Neighbor cell number
                dx = cell(nghbr_cell)%xc - cell(icell)%xc
                dy = cell(nghbr_cell)%yc - cell(icell)%yc
                dz = cell(nghbr_cell)%zc - cell(icell)%zc
                do mweight = 1,nlsq
                    lsq_weight_invdis_power = real(mweight-1, p2)
                    weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                    a(k,1,mweight) = weight_k*dx
                    a(k,2,mweight) = weight_k*dy
                    a(k,3,mweight) = weight_k*dz
                end do
                maxdx  = max(abs(dx),maxdx)
                maxdy  = max(abs(dy),maxdy)
                maxdz  = max(abs(dz),maxdz)
            end do nghbr_loop

            ! Do ghost cells if any
            nghbr_loop2 : do k = 1, lsqc(icell)%nbf
                ci = lsqc(icell)%gcells(1,k)
                ib = lsqc(icell)%gcells(2,k)
                dx = gcell(ib)%xc(ci) - cell(icell)%xc
                dy = gcell(ib)%yc(ci) - cell(icell)%yc
                dz = gcell(ib)%zc(ci) - cell(icell)%zc
                do mweight = 1,nlsq
                    lsq_weight_invdis_power = real(mweight-1, p2)
                    weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                    a(k+lsqc(icell)%n_nnghbrs,1,mweight) = weight_k*dx
                    a(k+lsqc(icell)%n_nnghbrs,2,mweight) = weight_k*dy
                    a(k+lsqc(icell)%n_nnghbrs,3,mweight) = weight_k*dz
                end do
                maxdx  = max(abs(dx),maxdx)
                maxdy  = max(abs(dy),maxdy)
                maxdz  = max(abs(dz),maxdz)
            end do nghbr_loop2
            !-------------------------------------------------------
            ! Perform QR factorization and compute R^{-1}*Q^T from A(m,n).
            do mweight = 1,nlsq
                call qr_factorization(a(:,:,mweight),rinvqt(:,:,mweight),m + lsqc(icell)%nbf,n)
            end do
            !-------------------------------------------------------
            ! Compute and store the LSQ coefficients: R^{-1}*Q^T*w.
            !
            ! (wx,wy,wz) = R^{-1}*Q^T*RHS
            !            = sum_k (cx,cy,cz)*(wk-wi).

            nghbr_loop3 : do k = 1, m
                nghbr_cell = lsqc(icell)%nghbr_lsq(k)
                dx = cell(nghbr_cell)%xc - cell(icell)%xc
                dy = cell(nghbr_cell)%yc - cell(icell)%yc
                dz = cell(nghbr_cell)%zc - cell(icell)%zc
                do mweight = 1,nlsq
                    lsq_weight_invdis_power = real(mweight-1, p2)
                    weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                    lsqc(icell)%cx(k,mweight)  = rinvqt(ix,k,mweight) * weight_k
                    lsqc(icell)%cy(k,mweight)  = rinvqt(iy,k,mweight) * weight_k
                    lsqc(icell)%cz(k,mweight)  = rinvqt(iz,k,mweight) * weight_k
                    lsqc(icell)%cf(:,k,mweight)  = rinvqt(:,k,mweight) * weight_k
                end do
            end do nghbr_loop3

            ! Do ghost cells if any
            nghbr_loop4 : do k = m+1, m + lsqc(icell)%nbf
                ci = lsqc(icell)%gcells(1,k - m)
                ib = lsqc(icell)%gcells(2,k - m)
                dx = gcell(ib)%xc(ci) - cell(icell)%xc
                dy = gcell(ib)%yc(ci) - cell(icell)%yc
                dz = gcell(ib)%zc(ci) - cell(icell)%zc
                do mweight = 1,nlsq
                    lsq_weight_invdis_power = real(mweight-1, p2)
                    weight_k = one / sqrt( dx**2 + dy**2 + dz**2 )**lsq_weight_invdis_power
                    lsqc(icell)%gcx(k - m,mweight)  = rinvqt(ix,k,mweight) * weight_k
                    lsqc(icell)%gcy(k - m,mweight)  = rinvqt(iy,k,mweight) * weight_k
                    lsqc(icell)%gcz(k - m,mweight)  = rinvqt(iz,k,mweight) * weight_k
                    lsqc(icell)%gcf(:,k - m,mweight)  = rinvqt(:,k,mweight) * weight_k
                end do
            end do nghbr_loop4
            !-------------------------------------------------------
            ! Deallocate a and rinvqt, whose size may change in the next cell. 
            deallocate(a, rinvqt)
        end do cloop

        ! Verification
        ! Compute the gradient of w = 2*x+y+4*z to se if we get wx = 2, wy = 1, and wz = 4 correctly
        verification_error = .false.

        do mweight = 1,nlsq
            do icell = 1,ncells
                ! initialize wx, wy, and wz
                wx = zero
                wy = zero
                wz = zero
                ! (xi,yi,zi) to be used to compute the function 2*x+y+4z at i
                xi = cell(icell)%xc
                yi = cell(icell)%yc
                zi = cell(icell)%zc

                ! look over the vertex neighboes
                do k = 1,lsqc(icell)%n_nnghbrs
                    nghbr_cell = lsqc(icell)%nghbr_lsq(k)
                    xk = cell(nghbr_cell)%xc
                    yk = cell(nghbr_cell)%yc
                    zk = cell(nghbr_cell)%zc
                    ! This is how we use the LSQ coefficients: accumulate cx*(wk-wi)
                    ! and cy*(wk-wi) and cz*(wk-wi)
                    wx = wx + lsqc(icell)%cx(k,mweight)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wy = wy + lsqc(icell)%cy(k,mweight)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wz = wz + lsqc(icell)%cz(k,mweight)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                end do
                do k = 1,lsqc(icell)%nbf
                    ci = lsqc(icell)%gcells(1,k)
                    ib = lsqc(icell)%gcells(2,k)
                    xk = gcell(ib)%xc(ci)
                    yk = gcell(ib)%yc(ci)
                    zk = gcell(ib)%zc(ci)
                    wx = wx + lsqc(icell)%gcx(k,mweight)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wy = wy + lsqc(icell)%gcy(k,mweight)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wz = wz + lsqc(icell)%gcz(k,mweight)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                end do

                maxDeltasNZ = zero
                if (maxdx > 0.001_p2) maxDeltasNZ(1) = one
                if (maxdy > 0.001_p2) maxDeltasNZ(2) = one
                if (maxdz > 0.001_p2) maxDeltasNZ(3) = one
                if ( maxDeltasNZ(1)*abs(wx-two) > 1.0e-06_p2 .or. &
                    maxDeltasNZ(2)*abs(wy-one) > 1.0e-06_p2 .or. &
                    maxDeltasNZ(3)*abs(wz-4.0_p2) > 1.0e-06_p2) then
                        write(*,*) " wx = ", wx, " exact ux = 2.0"!,maxDeltasNZ(1)*abs(wx-two)
                        write(*,*) " wy = ", wy, " exact uy = 1.0"!,maxDeltasNZ(2)*abs(wy-one)
                        write(*,*) " wz = ", wz, " exact uz = 4.0"!, maxDeltasNZ(3)*abs(wz-4.0_p2),maxDeltasNZ(3)
                        verification_error = .true.
                end if
            end do
        end do


        if (verification_error) then

            write(*,*) " LSQ coefficients are not correct. See above. Stop."
            stop
         
        else
         
            write(*,*) " Verified: LSQ coefficients are exact for a linear function."
         
        endif
         
        write(*,*)
        write(*,*) " End of Computing LSQ coefficients... "
        write(*,*) "--------------------------------------------------"
        write(*,*)
        
        ! End of Compute the LSQ coefficients in all cells.
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
         
    end subroutine compute_cell_coefficients

end module least_squares