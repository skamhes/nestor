module sort_routines
    ! A module containing various bespoke sorting routines to be used by other modules in this solver.

    use common , only : p2

    contains

    subroutine queued_natural_merge_sort(nn,nruns,runpointer,runlength,rtape,exclude,n_nnghbr,nnghbr)

        implicit none

        integer,               intent(in)    :: nn         ! length of rtape
        integer,               intent(in)    :: nruns      ! # of runs in the vector (also the number of attached nodes)
        integer, dimension(:), intent(inout) :: runpointer ! Pointer to the first # in each run
        integer, dimension(:), intent(inout) :: runlength  ! Length of each run
        integer, dimension(:), intent(inout) :: rtape      ! read tape
        integer,               intent(in)    :: exclude    ! number to be skipped in the merge (index of the host cell)

        integer,                        intent(out) :: n_nnghbr ! number of unique neighbors
        integer, dimension(:), pointer, intent(out) :: nnghbr

        ! Local vars
        integer, dimension(:), pointer  :: wtape        ! write tape
        integer                         :: irun, jrun   ! index of the two active runs
        integer                         :: h1, h2       ! read headers 1 and 2
        integer                         :: nr           ! nruns local
        integer                         :: lenwt        ! length of the write tape
        integer                         :: wh           ! write head
        integer                         :: ncloops

        allocate(wtape(nn))
        nr = nruns ! there could be unintended side effects if we change the var nruns points to
        wh = 0
        ncloops = 0
        

        cloop : do while(nruns > 1) 
            
            ncloops = ncloops + 1
            irun = 1
            jrun = 2

            gloop : do ! global loop not Augustus
                h1 = runpointer(irun)
                h2 = runpointer(jrun)
                rloop : do
                    ! Skip the exclude value
                    if (rtape(h1) == exclude) h1 = h1 + 1
                    if (rtape(h2) == exclude) h2 = h2 + 1

                    wh = wh + 1
                    ! Merge identical values
                    if (rtape(h1) == rtape(h2)) then
                        wtape(wh) = rtape(h1)
                        h1 = h1 + 1
                        h2 = h2 + 1
                    !Otherwise take the smaller value
                    elseif (rtape(h1) < rtape(h2)) then
                        wtape(wh) = rtape(h1)
                        h1 = h1 + 1
                    else ! rtape(h1) > rtape(h2)
                        wtape(wh) = rtape(h2)
                        h2 = h2 + 1
                    endif

                    ! Check if we are at the end of either run
                    if (h1 > runlength(irun)) then
                        do while(h2 <= runlength(jrun))
                            wh = wh + 1
                            wtape(wh) = rtape(h2)
                            h2 = h2 + 1
                        end do
                        exit rloop
                    end if
                    if (h2 > runlength(jrun)) then
                        do while(h1 <= runlength(irun))
                            wh = wh + 1
                            wtape(wh) = rtape(h1)
                            h1 = h1 + 1
                        end do
                        exit rloop
                    end if
                    
                end do rloop

                ! move the headers over 2
                irun = irun + 2
                jrun = jrun + 2

                if (irun > nruns) then
                     exit cloop ! all
                elseif (irun == nruns) then
                    h1 = runpointer(irun)
                    do while(h1 <= runlength(irun))
                        wh = wh + 1
                        wtape(wh) = rtape(h1)
                        h1 = h1 + 1
                    end do
                endif
            end do gloop
        
        end do cloop






    end subroutine queued_natural_merge_sort
end module sort_routines