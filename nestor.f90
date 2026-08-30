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
    use config, only : read_nml_config, generate_tec_file_b, write_restart, restart, line_implicit

    use common, only : version

    use files,  only : set_filenames
    
    use inout,  only : write_tecplot_file_b, write_restart_file

    use grid,   only : read_grid, read_su2, construct_grid, cell, ncells, face, nfaces, nb, bound

    use reorder , only : reorder_rcm

    use limplicit , only : build_lines

    use solution, only : allocate_solution_vars, define_problem

    use steady_solver, only : steady_solve

    use utils , only : isolver_type, SOLVER_IMPLICIT, isolver_type, SOLVER_GCR

    use initialize , only : init_jacobian, set_initial_solution

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

    call define_problem

    call construct_grid

    call reorder_rcm

    if (line_implicit) call build_lines

    if (isolver_type == SOLVER_IMPLICIT .OR. isolver_type == SOLVER_GCR ) call init_jacobian

    call allocate_solution_vars

    call set_initial_solution

    call steady_solve

    if (write_restart) then
        call write_restart_file
    end if

    if ( generate_tec_file_b ) then
        call write_tecplot_file_b
    end if
    
    ! if ( generate_tec_file_v ) then
    !     call write_tecplot_file_v
    ! end if
end program nestor