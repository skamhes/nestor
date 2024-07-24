module files

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
end module files