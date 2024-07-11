module gradient

    use common , only : p2

    implicit none 

    public init_gradients

    contains

    subroutine init_gradients
        
        use least_squares , only : construct_lsq_stencil

        use config        , only : grad_method, lsq_stencil

        use solution      , only : ccgradq, vgradq

        use common        , only : p2, zero

        !use grid          , only : nnodes, ncells


        implicit none
        
        if (trim(grad_method) == 'lsq') then
            ! Build LSQ Stencil and coefficients
            call construct_lsq_stencil
            if (trim(lsq_stencil) == 'w_vertex') then
                ! Gradient arrays are allocated in solution subroutine
                write(*,*)
                write(*,*) 'Initializing gradient vertex arrays.'
                vgradq = zero
            endif
        endif

        ! Initialize the cell centered array always
        write(*,*)
        write(*,*) 'Initializing gradient cell centered arrays.'
        ccgradq = zero

    end subroutine init_gradients
end module gradient