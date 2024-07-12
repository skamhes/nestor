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

        use gradient        , only : com

        implicit none

        ! Grid Vars
        real(p2)                    :: xm, ym, zm
        integer                     :: c1, c2,  v1, v2, v3,
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2

        ! Flow variables
        real(p2), dimension(5)      :: u1, u2, q1, q2
        real(p2), dimension(3,5)    :: gradq1, gradq2, gradqb
        real(p2), dimension(5)      :: num_flux
        real(p2), dimension(5)      :: qb
        real(p2)                    :: wave_speed
        real(p2)                    :: phi1, phi2

        ! Misc int/counters
        integer                     :: i, os
        integer                     :: j, ib, ix, iu, ii
        
        !--------------------------------------------------------------------------------
        ! Initialize the residuals and wsn = the sum of (max_wave_speed)*(face length)).
        ! Note: wsn is required to define a time step.
        cell_loop1 :  do i = 1, ncells

            res(:,i) = zero
            wsn(i)   = zero
     
        end do cell_loop1

        !--------------------------------------------------------------------------------
        ! Compute gradients at cells.
        !
        if (accuracy_order == 2) then
            call 

    end subroutine compute_residual

end module residual