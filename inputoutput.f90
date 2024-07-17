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
                                     bound

        use solution        , only : q, nq

        use config          , only : project_name
                                     
        implicit none 

        integer :: i, os, ibn
        
        integer                           :: j, k, ib, bcell_i, candidate_node, nj, ni
        real(p2), dimension(:,:), pointer :: qn
        integer , dimension(:  ), pointer :: nc
        logical, dimension(:), pointer    :: already_added
        real(p2)                          :: an
        integer, dimension(:)             :: nbnodes

        allocate(qn(nq, nnodes))
        allocate(nc(    nnodes))
        
        allocate(already_added(nnodes))
        allocate(nbnodes(nb))

        
        nc = 0
        qn = zero

        bound_loop : do ib = 1,nb
            bface_loop : do i = 1,bound(ib)%nbfaces
                bcell_i = bound(ib)% bcell(i)
                do k = 1,cell(bcell_i)%nvtx
                    qn(:,cell(bcell_i)%vtx(k)) = qn(:,cell(bcell_i)%vtx(k)) + q (:,bcell_i)
                    nc(  cell(bcell_i)%vtx(k)) = nc(  cell(bcell_i)%vtx(k)) + 1
                enddo
            enddo bface_loop
        enddo bound_loop
        do i = 1,nnodes
            qn(:,i) = qn(:,i) / nc(i)
            Mn(j) = sqrt(wn(2,j)**2 + wn(3,j)**2 + wn(4,j)**2)  ! mach number (wrt free stream a)
            an    = sqrt(gamma*wn(5,j)/wn(1,j))                   ! local speed of sound
            Mn(j) = Mn(j) / an  
        end do

        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) ' Writing Tecplot file = ', trim(filename_tecplot_b)
        write(*,*)
    
        !Open the output file.
        open(unit=8, file=filename_tecplot_b, status="unknown", iostat=os)   

        !---------------------------------------------------------------------------

        !(0)Header information

        write(8,*) 'TITLE = "GRID"'
        write(8,*) 'VARIABLES = "x","y","z","rho","u","v","w","p","M"'

        do ib = 1,nb
            already_added = .false.
            do j = 1,bound(ib)%nbfaces
                nj = bound(ib)%bfaces(1,j)
                do i = 2,nj
                    ni = bound(ib)%bfaces(i,j)
                    if (.not.already_added(ni)) then
                        already_added(ni) = .true.
                        nbnodes(ib) = nbnodes(ib) + 1
                    endif
                enddo
            enddo
        enddo
        do ib = 1,nb
            already_added = .false.
            write(8,*) 'ZONE T = "',ib,'"  n=', nbnodes(ib), &
                        ',e=', bound(ib)%nbfaces,' , zonetype=fequadrilateral, datapacking=point'

            ! Loop through nodes
            do j = 1,bound(ib)%nbfaces
                nj = bound(ib)%bfaces(1,j)
                do i = 2,nj
                    ni = bound(ib)%bfaces(i,j)
                    if (.not.already_added(ni)) then
                        already_added(ni) = .true.
                        write(8,'(11es25.15)') x(ni), y(ni), z(ni), wn(ir,ni), wn(iu,ni), wn(iv,ni), wn(iw,ni), wn(ip,ni), Mn(ni)
                    endif
                enddo
            enddo

            ! Loop through faces
            do i = 1,bound(ib)%nbfaces
                if (bound(ib)%bfaces(1,i) == 3) then ! write tri as a degenerate quad
                    write(8,'(4i10)') bound(ib)%bfaces(2,i), bound(ib)%bfaces(3,i), & 
                                        bound(ib)%bfaces(4,i), bound(ib)%bfaces(4,i)
                else if (bound(ib)%bfaces(1,i) == 4) then ! write quad
                    write(8,'(4i10)') bound(ib)%bfaces(2,i), bound(ib)%bfaces(3,i), & 
                    bound(ib)%bfaces(4,i), bound(ib)%bfaces(5,i)
                end if
            enddo
        end do
        close(8)
        write(*,*)
        write(*,*) ' End of Writing Tecplot file = ', trim(filename_tecplot_b)
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    end subroutine write_tecplot_file_b
end module inout