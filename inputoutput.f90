module inout

    implicit none

    public

    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    ! File names

    character(80) :: filename_grid         ! input grid filename (.ugrid)
    character(80) :: filename_bc           ! input bc   filename (.bcmap)
    character(80) :: filename_tecplot_b    ! output tecplot boundary filename (.dat)
    character(80) :: filename_tecplot_v    ! output tecplot volume filename (.dat)
    character(80) :: filename_data         ! Output of U array
    
    public :: write_tecplot_file_b
    ! public :: write_tecplot_file_v ! TODO

    contains 

    subroutine set_filenames
        use config , only : project_name, grid_type

        implicit none
        
        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) " Setting up file names..... "
        write(*,*)
    
        !-----------------------------------------------------------------------
        ! Input grid file (.ugrid):
        ! E.g., filename_grid = "test.grid" if project_name = "test".
        if (trim(grid_type) == 'ugrid') then
            filename_grid   = trim(project_name) // '.ugrid'
        elseif (trim(grid_type) == 'su2') then
            filename_grid   = trim(project_name) // '.su2'   
        else
            write(*,*) 'Unsupported grid type: ', trim(grid_type), '. Stop!'
            stop
        endif
        write(*,'(a28,a28)') "   filename_grid = ", trim(filename_grid)
    
        !-----------------------------------------------------------------------
        ! Input boundary condition file (ASCII file)
        ! E.g., filename_bc = "test.bc" if project_name = "test".
    
            filename_bc   = trim(project_name) // '.bcmap'
    
            write(*,'(a28,a28)') "   filename_bc = ", trim(filename_bc)
    
        !-----------------------------------------------------------------------
        ! Output: Tecplot boundary file (ASCII file)
        ! E.g., filename_tecplot = "test_tec.dat" if project_name = "test".
    
            filename_tecplot_b = trim(project_name) // '_b_tec.dat'
    
            write(*,'(a28,a28)') " filename_tecplot_b = ", trim(filename_tecplot_b)

        !-----------------------------------------------------------------------
        ! Output: Tecplot volume file (ASCII file)
        ! E.g., filename_tecplot = "test_tec.dat" if project_name = "test".
    
            filename_tecplot_v = trim(project_name) // '_v_tec.dat'
    
            write(*,'(a28,a28)') " filename_tecplot_b = ", trim(filename_tecplot_b)


        !-----------------------------------------------------------------------
        ! Output: Tecplot boundary file (ASCII file)
        ! E.g., filename_tecplot = "test_tec.dat" if project_name = "test".
    
            filename_data = trim(project_name) // '.kdat'
    
            write(*,'(a28,a28)') "       filename_data = ", trim(filename_data)

    
        write(*,*)
        write(*,*) " End of Setting up file names..... "
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    
    end subroutine set_filenames
    
    subroutine write_tecplot_file_b
        use common          , only : p2, zero
                                     
        use grid            , only : nnodes, x, y, z, &
                                     ntria, tria, &
                                     nquad, quad, &
                                     bc_type, nb, &
                                     cell, ncells, &
                                     bound, bgrid_type

        use solution        , only : q, nq, gamma

        use config          , only : project_name

        use lowlevel        , only : my_alloc_int_ptr
                                     
        implicit none 

        type bnode_type
                integer                               :: nbnodes 
                integer, dimension(:),  pointer       :: bnodes  
        end type bnode_type
        type(bnode_type), dimension(nb)               :: bnode_data
            
        type(bgrid_type), dimension(:), pointer       :: bound_export
        integer :: i, os, ibn
        
        integer                           :: j, k, ib, nk, nj, ni, j_count
        integer                           :: bcell_i, candidate_node
        real(p2), dimension(:,:), pointer :: qn
        real(p2), dimension(:), pointer   :: Mn, rhon
        integer , dimension(:  ), pointer :: nc
        logical                           :: already_added
        real(p2)                          :: an
        ! integer, dimension(:), pointer    :: nbnodes

        allocate(qn(nq, nnodes))
        allocate(nc(    nnodes))
        
        allocate(Mn(           nnodes))
        allocate(rhon(         nnodes))

        
        nc = 0
        qn = zero

        bound_loop : do ib = 1,nb
            bface_loop : do i = 1,bound(ib)%nbfaces
                bcell_i = bound(ib)%bcell(i)
                do k = 1,cell(bcell_i)%nvtx
                    qn(:,cell(bcell_i)%vtx(k)) = qn(:,cell(bcell_i)%vtx(k))  + q(:,bcell_i)
                    nc(cell(bcell_i)%vtx(k))   = nc(cell(bcell_i)%vtx(k)) + 1
                end do
            end do bface_loop
        end do bound_loop

        do j = 1,nnodes
            qn(:,j) = qn(:,j) / nc(j) ! copmute an average
            Mn(j) = sqrt(qn(2,j)**2 + qn(3,j)**2 + qn(4,j)**2)  ! mach number (wrt free stream a)
            an    = sqrt(gamma*qn(5,j)/qn(1,j))                 ! local speed of sound
            Mn(j) = Mn(j) / an                                  ! local mach number
            ! rho    = p      * gamma / T
            rhon(j) = qn(1,j) * gamma / qn(5,j)
        end do

        allocate(bound_export(nb))
        boundary_loop : do ib = 1,nb
            bnode_data(ib)%nbnodes = zero
            allocate(bnode_data(ib)%bnodes(1))
            bnode_data(ib)%bnodes(1) = zero
            bound_export(ib)%nbfaces = bound(ib)%nbfaces
            allocate( bound_export(ib)%bfaces( 5,bound(ib)%nbfaces ) )  
            bound_export(ib)%bfaces = zero
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                bound_export(ib)%bfaces(1,i) = bound(ib)%bfaces(1,i)
                bface_vertex_loop : do k = 2,(bound(ib)%bfaces(1,i) + 1) ! loop through number of vertices for face
                    if (bnode_data(ib)%nbnodes == 0 ) then
                        bnode_data(ib)%nbnodes = bnode_data(ib)%nbnodes + 1
                        bnode_data(ib)%bnodes  = bound(ib)%bfaces(k,i)
                        bound_export(ib)%bfaces(k,i) = 1
                        cycle bface_vertex_loop 
                    end if 
                    candidate_node = bound(ib)%bfaces(k,i)
                    already_added = .false.
                    bnodes_loop : do nk = 1,bnode_data(ib)%nbnodes
                        if (candidate_node == bnode_data(ib)%bnodes(nk)) then
                            already_added = .true.
                            bound_export(ib)%bfaces(k,i) = nk
                            exit bnodes_loop  
                        end if 
                    end do bnodes_loop 
                    if (.not.already_added) then
                        bnode_data(ib)%nbnodes = bnode_data(ib)%nbnodes + 1
                        call my_alloc_int_ptr(bnode_data(ib)%bnodes,bnode_data(ib)%nbnodes)
                        bnode_data(ib)%bnodes(bnode_data(ib)%nbnodes) = candidate_node
                        bound_export(ib)%bfaces(k,i) = bnode_data(ib)%nbnodes
                    end if
                end do bface_vertex_loop 
            end do bfaces_loop
        end do boundary_loop

        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) ' Writing Tecplot file = ', trim(filename_tecplot_b)
        write(*,*)
    
        !Open the output file.
        open(unit=8, file=filename_tecplot_b, status="unknown", iostat=os)   

        !---------------------------------------------------------------------------

        !(0)Header information

        write(8,*) 'TITLE = "GRID"'
        write(8,*) 'VARIABLES = "x","y","z","p","u","v","w","T","rho","M"'

        do ib = 1,nb
            write(8,*) 'ZONE T = "',ib,'"  n=', bnode_data(ib)%nbnodes, &
                            ',e=', bound(ib)%nbfaces,' , zonetype=fequadrilateral, datapacking=point'
            do j_count = 1,bnode_data(ib)%nbnodes
                j = bnode_data(ib)%bnodes(j_count)
                write(8,'(10es25.15)') x(j), y(j), z(j), qn(1,j), qn(2,j), qn(3,j), qn(4,j), qn(5,j), rhon(j), Mn(j)
            end do
            ! Loop through faces
            do i = 1,bound_export(ib)%nbfaces
                if (bound_export(ib)%bfaces(1,i) == 3) then ! write tri as a degenerate quad
                    write(8,'(4i10)') bound_export(ib)%bfaces(2,i), bound_export(ib)%bfaces(3,i), & 
                                        bound_export(ib)%bfaces(4,i), bound_export(ib)%bfaces(4,i)
                else if (bound_export(ib)%bfaces(1,i) == 4) then ! write quad
                    write(8,'(4i10)') bound_export(ib)%bfaces(2,i), bound_export(ib)%bfaces(3,i), & 
                    bound_export(ib)%bfaces(4,i), bound_export(ib)%bfaces(5,i)
                end if

            end do
        end do
        
        close(8)
        write(*,*)
        write(*,*) ' End of Writing Tecplot file = ', trim(filename_tecplot_b)
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    end subroutine write_tecplot_file_b
end module inout