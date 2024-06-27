!********************************************************************************
! NaviEr-STOkes Resolved.  Why that name?  I wanted to name a code after my cat...
! He gets something names after him, I get someone to blame for any bugs (I swear 
! he walked across the keyboard...)
!
! Author: Karsten Hendrickson
! Version: 0.0.1

program nestor
    use config, only : read_nml_config, grid_type

    use common, only : version

    use inout,  only : set_filenames

    use grid,   only : read_grid, read_su2, construct_grid

    use solution, only : allocate_solution_vars

    use steady_solver, only : steady_solve

    implicit none

    write(*,*)
    write(*,*) "----------------------------------------------------------------"
    write(*,*)
    write(*,'(a,3(i1,a))') "  Nestor Version: ", version(1),".",version(2),".",version(3),"."
    write(*,*)
    write(*,*) "----------------------------------------------------------------"
    write(*,*)   

    !-------------------------------------------------------------------------------
    ! READ CONFIG SETTINGS
    !-------------------------------------------------------------------------------
    call read_nml_config("nestor.nml")

    !-------------------------------------------------------------------------------
    ! DEFINE INPUT AND OUTPUT FILENAMES
    !-------------------------------------------------------------------------------
    call set_filenames

    !-------------------------------------------------------------------------------
    ! READ GRID
    !-------------------------------------------------------------------------------
    if (trim(grid_type) == 'ugrid') then
        call read_grid
    elseif (trim(grid_type) == 'su2') then
        call read_su2
    else
        write(*,*) 'Unsupported grid type: ', trim(grid_type), '. Stop!'
        stop
    endif

    call construct_grid

    call allocate_solution_vars

    call steady_solve
end program nestor