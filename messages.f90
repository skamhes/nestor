module messages

    ! this module handles sending messages to the console.  This is anything that wouldn't be output in a "normal" or "well behaved"
    ! run.  (therefore most of goes to stderr)

    use iso_fortran_env

    implicit none

    public

    contains


    subroutine print_error_location(name_subroutine, line_number, filename, stream)

        implicit none

        character(len=*), intent(in) :: name_subroutine, filename
        integer,          intent(in) :: line_number
        integer, optional,intent(in) :: stream ! file descriptor

        integer :: fd = ERROR_UNIT ! Default file descriptor is stderr
        character(80) :: cline_number

        write(cline_number,'(I0)') line_number

        if (present(stream)) fd = stream

        write(*,*) ' Error occured in ',trim(name_subroutine),' in ',trim(filename),':',trim(cline_number),'. Stopping...'
    end subroutine print_error_location
    
end module messages