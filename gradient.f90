module gradient

    use common , only : p2

    implicit none 

    public init_gradients

    contains

    subroutine init_gradients
        
        use least_squares , only : construct_lsq_stencil

        use config        , only : grad_method

        implicit none
        
        call construct_lsq_stencil

        ! call compute_lsq_weights
    end subroutine init_gradients
end module gradient