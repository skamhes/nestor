!********************************************************************************
! NaviEr-STOkes Resolved.  Why that name?  I wanted to name a code after my cat...
! He gets something names after him, I get someone to blame for any bugs (I swear 
! he walked across the keyboard...)
!
! Author: Karsten Hendrickson
! Version: 0.0.1

program nestor
    use config, only : read_nml_config

    use common, only : version

    use inout,  only : set_filenames
    implicit none

    write(*,*)
    write(*,*) "----------------------------------------------------------------"
    write(*,*)
    write(*,*) "  Nestor Version: ", version(1),".",version(2),".",version(3),"."
    write(*,*)
    write(*,*) "----------------------------------------------------------------"
    write(*,*)   


    call read_nml_config("nestor.nml")

    set_filenames
end program nestor