module residual

    implicit none 
    
    public :: compute_residual

    contains 

    subroutine compute_residual

        use common          , only : p2, zero, half, one, two

        use config          , only : method_inv_flux, accuracy_order, use_limiter

        use grid            , only : ncells, cell,  &
                                     nfaces, face,  &
                                     nb,     bound, &
                                     bc_type,       &
                                     face_nrml,     &
                                     face_nrml_mag, &
                                     face_centroid
        
        use solution        , only : res, q, ccgradq, wsn, q2u

        use inviscid_flux   , only : compute_inviscid_flux

        use limiter         , only : compute_limiter

        use bc_states       , only : get_right_state

    end subroutine compute_residual

end module residual