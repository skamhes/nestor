!********************************************************************************
! NaviEr-STOkes Robust methods.  Why that name?  I wanted to name a code after my cat...
! He gets something names after him, I get someone to blame for any bugs (I swear 
! he walked across the keyboard...)  Also, once I have a working code I think it would 
! be an interesting project to focus on implementing numerical methods that focus on robust 
! solutions...
! 
!
! Author: Karsten Hendrickson
! Version: 0.0.1

program nestor
    use config, only : read_nml_config, generate_tec_file_b

    use common, only : version

    use files,  only : set_filenames
    
    use inout,  only : write_tecplot_file_b

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
    call read_grid

    call construct_grid

    call allocate_solution_vars

    call steady_solve

    ! if (write_data) then
    !     call write_data_file
    ! end if

    if ( generate_tec_file_b ) then
        call write_tecplot_file_b
    end if
    
    ! if ( generate_tec_file_v ) then
    !     call write_tecplot_file_v
    ! end if
end program nestor