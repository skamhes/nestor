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

    type lsq_data_type
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
    end type lsq_data_type

    !Cell data array in the custom data type.
    type(lsq_data_type), dimension(:), pointer :: lsq  !cell-centered LSQ array

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

        use utils , only : ilsq_stencil, LSQ_STENCIL_WVERTEX

        implicit none

        select case(ilsq_stencil)
        case(LSQ_STENCIL_WVERTEX)
            call construct_wvertex_stencil
            call compute_vertex_coefficients
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
        allocate( lsq(nnodes ) )

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
            lsq(i)%ncells_lsq = node(i)%nic + node(i)%nbc
            allocate( lsq(i)%cell_lsq( lsq(i)%ncells_lsq  ) )
            allocate( lsq(i)%ib_lsq( lsq(i)%ncells_lsq ) )
            lsq(i)%btype = 0 ! stays zero for internal cells

            node(i)%nic = 0
            node(i)%nbc = 0
        end do

        ! Add internal cells to the array
        do i = 1,ncells
            do k = 1, cell(i)%nvtx
                vk           = cell(i)%vtx(k)
                node(vk)%nic = node(vk)%nic + 1
                lsq(vk)%cell_lsq(node(vk)%nic) = i
                lsq(vk)%ib_lsq(node(vk)%nic) = INTERNAL ! Initialize the cell as internal
            end do
        end do

        ! Add boundary cells to the array
        do ib = 1,nb
            do j = 1,bound(ib)%nbfaces
                fj           = bound(ib)%bfaces(1,j) ! (/#_of_nodes v1:v_end/)
                do k = 2,(fj+1) 
                    vk           = bound(ib)%bfaces(k,j)
                    node(vk)%nbc = node(vk)%nbc + 1
                    lsq(vk)%cell_lsq(node(vk)%nic + node(vk)%nbc) = j
                    ! We will need face data in some cases so we will use bface# in place of cell
                    ! In the future I can sort these so the boundary faces match with their attached cell
                    ! In theory this would allow me to skip some assignmets.  Will do another time.
                    
                    ! Give boundary of attached bcell
                    lsq(vk)%ib_lsq(node(vk)%nic + node(vk)%nbc) = ib ! internal cell

                    ! Assign the boundary condition based on above hierarchy
                    select case(trim(bc_type(ib))) 
                        case('freestream')
                            lsq(vk)%btype = max(lsq(vk)%btype,FREE_STREAM)
                        case('symmetry')
                            ! symmetry is treated numerically like a slip_wall
                            lsq(vk)%btype = max(lsq(vk)%btype,SLIP_WALL)
                        case('slip_wall')
                            lsq(vk)%btype = max(lsq(vk)%btype,SLIP_WALL)
                        case('no_slip_wall')
                            lsq(vk)%btype = max(lsq(vk)%btype,NO_SLIP_WALL)
                        case('outflow_subsonic')
                            lsq(vk)%btype = max(lsq(vk)%btype,PRESSURE_OUTLET)
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
            allocate( lsq(i)%cx4(lsq(i)%ncells_lsq) )
            allocate( lsq(i)%cy4(lsq(i)%ncells_lsq) )
            allocate( lsq(i)%cz4(lsq(i)%ncells_lsq) )
            allocate( lsq(i)%cq4(lsq(i)%ncells_lsq) )
            if ( lsq(i)%btype > 0 ) then
                allocate( lsq(i)%cx3(lsq(i)%ncells_lsq) )
                allocate( lsq(i)%cy3(lsq(i)%ncells_lsq) )
                allocate( lsq(i)%cz3(lsq(i)%ncells_lsq) )
            endif
        end do
        ! I want the option to prescribe the node value at boundaries where it's known.
        ! It seems like this should provide an increased robustness(?).  Intuition says I can't
        ! just use the calculated weights for cy, cx, and cz, but I'll have to investigate...
        
    end subroutine construct_wvertex_stencil

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
            m = lsq(i)%ncells_lsq ! # of connected cells
            allocate( a3(m,3) )
            allocate( a4(m,4) )
            allocate( rinvqt3(3,m) )
            allocate( rinvqt4(4,m) )

            connect_loop : do k = 1,m
                if ( lsq(i)%ib_lsq(k) == INTERNAL ) then ! Internal ib = 0
                    connect_cell = lsq(i)%cell_lsq(k)
                    dx = cell(connect_cell)%xc - x(i)
                    dy = cell(connect_cell)%yc - y(i)
                    dz = cell(connect_cell)%zc - z(i)
                else
                    connect_bface = lsq(i)%cell_lsq(k)
                    ib            = lsq(i)%ib_lsq(k)
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
                if ( lsq(i)%ib_lsq(k) == INTERNAL ) then ! Internal ib = 0
                    connect_cell = lsq(i)%cell_lsq(k)
                    dx = cell(connect_cell)%xc - x(i)
                    dy = cell(connect_cell)%yc - y(i)
                    dz = cell(connect_cell)%zc - z(i)
                else
                    connect_bface = lsq(i)%cell_lsq(k)
                    ib            = lsq(i)%ib_lsq(k)
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
                lsq(i)%cx3(k) = rinvqt3(1,k) * weight_k
                lsq(i)%cy3(k) = rinvqt3(2,k) * weight_k
                lsq(i)%cz3(k) = rinvqt3(3,k) * weight_k
                ! 4 unknowns
                lsq(i)%cq4(k) = rinvqt4(1,k) * weight_k
                lsq(i)%cx4(k) = rinvqt4(2,k) * weight_k
                lsq(i)%cy4(k) = rinvqt4(3,k) * weight_k
                lsq(i)%cz4(k) = rinvqt4(4,k) * weight_k
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
            do k = 1,lsq(i)%ncells_lsq
                if ( lsq(i)%ib_lsq(k) == INTERNAL ) then
                    connect_cell = lsq(i)%cell_lsq(k)
                    xk = cell(connect_cell)%xc
                    yk = cell(connect_cell)%yc
                    zk = cell(connect_cell)%zc
                else
                    connect_bface = lsq(i)%cell_lsq(k) 
                    ib            = lsq(i)%ib_lsq(k)
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
                    wx = wx + lsq(i)%cx4(k)*( (2.0*xk+yk+4.0*zk) )
                    wy = wy + lsq(i)%cy4(k)*( (2.0*xk+yk+4.0*zk) )
                    wz = wz + lsq(i)%cz4(k)*( (2.0*xk+yk+4.0*zk) )
                    wq = wq + lsq(i)%cq4(k)*( (2.0*xk+yk+4.0*zk) )
                    ! we don't need q
                else ! unknowns_ == 3
                    wx = wx + lsq(i)%cx3(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yk+4.0*zi) )
                    wy = wy + lsq(i)%cy3(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                    wz = wz + lsq(i)%cz3(k)*( (2.0*xk+yk+4.0*zk)-(2.0*xi+yi+4.0*zi) )
                endif
            end do
            ! Loop through attached cells...again
            cell_grad_loop : do k = 1,lsq(i)%ncells_lsq
                if ( lsq(i)%ib_lsq(k) == INTERNAL ) then
                    connect_cell = lsq(i)%cell_lsq(k)
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

end module least_squares