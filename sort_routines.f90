module sort_routines
    ! A module containing various bespoke sorting routines to be used by other modules in this solver.

    use common , only : p2

    contains

    subroutine queued_natural_merge_sort(nn,nruns,runpointer,rtape,exclude,n_nnghbr,nnghbr)

        use lowlevel , only : ip_swap
        implicit none

        integer,                        intent(in)    :: nn         ! length of rtape
        integer,                        intent(in)    :: nruns      ! # of runs in the vector (also the number of attached nodes)
        integer, dimension(:), pointer, intent(inout) :: runpointer ! Pointer to the first # in each run
        integer, dimension(:), pointer, intent(inout) :: rtape      ! read tape
        integer,                        intent(in)    :: exclude    ! number to be skipped in the merge (index of the host cell)

        integer,                        intent(out)   :: n_nnghbr ! number of unique neighbors
        integer, dimension(:), pointer, intent(out)   :: nnghbr

        ! Local vars
        integer, dimension(:), pointer  :: wtape        ! write tape
        integer, dimension(:), pointer  :: wpointer     ! write tape start index for each written run
        integer                         :: irun, jrun   ! index of the two active runs
        integer                         :: h1, h2       ! read headers 1 and 2
        integer                         :: nr           ! nruns local
        integer                         :: lenwt        ! length of the write tape
        integer                         :: wh           ! write head
        integer                         :: ngloops

        allocate(wtape(nn))
        allocate(wpointer(nruns + 1))
        wpointer(1) = 1
        nr = nruns ! there could be unintended side effects if we change the var nruns points to
        ngloops = 0
        

        ! Not super thrilled about the amount of nesting
        cloop : do  ! should run ceil(nruns/2) times
            wh = 0
            
            irun = 1
            jrun = 2

            gloop : do ! global loop not Augustus. Should run nruns times
                ngloops = ngloops + 1
                h1 = runpointer(irun)
                h2 = runpointer(jrun)
                rloop : do ! # of loops equal 
                    ! Check if we are at the end of either run
                    if (h1 >= runpointer(irun+1)) then
                        ! call copy_remaining(h2,wh,runpointer(jrun+1),exclude,rtape,wtape)
                        do while(h2 < runpointer(jrun+1))
                            if (rtape(h2) /= exclude) then
                                wh = wh + 1
                                wtape(wh) = rtape(h2)
                            endif
                            h2 = h2 + 1
                        end do
                        exit rloop
                    end if
                    if (h2 >= runpointer(jrun+1)) then
                        do while(h1 < runpointer(irun+1))
                            if (rtape(h1) /= exclude) then
                                wh = wh + 1
                                wtape(wh) = rtape(h1)
                            endif
                            h1 = h1 + 1
                        end do
                        exit rloop
                    end if

                    ! Skip the exclude value
                    if (rtape(h1) == exclude) then 
                        h1 = h1 + 1
                        if (h1 >= runpointer(irun+1)) then
                            lh1 : do while(h2 < runpointer(jrun+1))
                                if (rtape(h2) /= exclude) then
                                    wh = wh + 1
                                    wtape(wh) = rtape(h2)
                                endif
                                h2 = h2 + 1
                            end do lh1
                            exit rloop
                        end if
                    endif
                    if (rtape(h2) == exclude) then 
                        h2 = h2 + 1
                        if (h2 >= runpointer(jrun+1)) then
                            lh2 : do while(h1 < runpointer(irun+1))
                                if (rtape(h1) /= exclude) then
                                    wh = wh + 1
                                    wtape(wh) = rtape(h1)
                                end if
                                h1 = h1 + 1
                            end do lh2
                            exit rloop
                        end if
                    endif

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
                end do rloop

                wpointer(ngloops + 1) = wh + 1

                ! move the headers over 2
                irun = irun + 2
                jrun = jrun + 2

                if (irun > nr) then
                     exit gloop ! all
                elseif (irun == nr) then
                    h1 = runpointer(irun)
                    do while(h1 < runpointer(irun+1))
                        if (rtape(h1) /= exclude) then
                            wh = wh + 1
                            wtape(wh) = rtape(h1)
                        end if
                        h1 = h1 + 1
                    end do
                    wpointer(ngloops + 1) = wh + 1
                    exit gloop
                endif
            end do gloop
        
            nr = ngloops

            if (nr < 2) exit cloop

            call ip_swap(wpointer,runpointer)
            call ip_swap(wtape,   rtape)

            wpointer = 0 ! not necessary since the length of data to be written shrinks with each iteration
            wtape    = 0 ! but helpful for debugging
            wpointer(1) = 1
            ngloops = 0

        end do cloop

        n_nnghbr = wh
        allocate(nnghbr(n_nnghbr))
        nnghbr(:) = wtape(1:n_nnghbr)




    end subroutine queued_natural_merge_sort

    subroutine copy_remaining(rh,wh,rpp1,exc,rt,wt)

        implicit none

        integer,                        intent(inout) :: rh     ! read header
        integer,                        intent(inout) :: wh     ! write header
        integer,                        intent(in)    :: rpp1   ! run pointer(i+1). Pointer to the start of the next run
        integer,                        intent(in)    :: exc    ! exclude value
        integer, dimension(:), pointer, intent(in)    :: rt     ! read tape
        integer, dimension(:), pointer, intent(inout) :: wt     ! write tape

        do while(rh < rpp1)
            if (rt(rh) /= exc) then
                wh = wh + 1
                wt(wh) = rt(rh)
            endif
            rh = rh + 1
        end do
    end subroutine copy_remaining
end module sort_routines