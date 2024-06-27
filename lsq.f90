module least_squares
    
    use common , only : p2

    implicit none

    type lsq_data_type
        integer                           ::   ncells_lsq  ! number of cells attached to lsq vertex
        integer , dimension(:)  , pointer ::     cell_lsq  ! lsq cell
        integer , dimension(:)  , pointer ::       ib_lsq  ! boundary number for each attached cell (0 = internal cell)
        real(p2), dimension(:)  , pointer ::           cq  ! LSQ coefficient for q at vertex
        real(p2), dimension(:)  , pointer ::           cx  ! LSQ coefficient for x-derivative
        real(p2), dimension(:)  , pointer ::           cy  ! LSQ coefficient for y-derivative
        real(p2), dimension(:)  , pointer ::           cZ  ! LSQ coefficient for z-derivative
        integer                           ::        btype  ! See below for values. 0 = internal cell
    end type lsq_data_type

    !Cell data array in the custom data type.
    type(lsq_data_type), dimension(:), pointer :: lsq  !cell-centered LSQ array

    ! Boundary int rankings
    integer, parameter :: NO_SLIP_WALL = 99
    integer, parameter :: SLIP_WALL = 98
    integer, parameter :: FREE_STREAM = 50
    integer, parameter :: PRESSURE_OUTLET = 20

    public construct_lsq_stencil

    contains

    subroutine construct_lsq_stencil

        use common , only : p2

        use config , only : lsq_stencil

        implicit none

        if (trim(lsq_stencil) == "w_vertex") then
            call construct_wvertex
        else
            write(*,*) "Unsupported LSQ Stencil"
            stop
        endif



    end subroutine construct_lsq_stencil

    subroutine construct_wvertex

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
        integer :: ck, vk, fj

        write(*,*)
        write(*,*) " --------------------------------------------------"
        write(*,*) " Constructing vertex neighbors... "
        write(*,*)

        allocate( node(nnodes) )

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
            lsq(i)%ncells_lsq = node(i)%nic + node(i)%nic
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
                lsq(vk)%ib_lsq(node(vk)%nic) = 0 ! internal cell
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
        
    end subroutine construct_wvertex


end module least_squares