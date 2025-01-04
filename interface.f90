module interface

    implicit none
    
    contains
    
    subroutine interface_flux(q1, q2, gradq1, gradq2, n12, &
                             xc1, yc1, zc1, xc2, yc2, zc2, &
                     xm, ym, zm, phi1, phi2, num_flux, wsn )

        use common                 , only : p2

        use config                 , only : accuracy_order

        use utils                  , only : imethod_inv_flux, IFLUX_ROE

        use solution               , only : q2u

        use inviscid_flux          , only : roe
        implicit none

        ! Inputs
        real(p2), dimension(5),     intent(in) :: q1, q2            ! Primative vars (prim for low mach)
        real(p2), dimension(3,5),   intent(in) :: gradq1, gradq2    ! Gradients of primitive vars
        real(p2), dimension(3),     intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                   intent(in) :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                   intent(in) :: xc2, yc2, zc2     ! Right cell centroid
        real(p2),                   intent(in) :: xm, ym, zm        ! Face midpoint
        real(p2),                   intent(in) :: phi1, phi2        ! Limiter

        ! Output
        real(p2), dimension(5),     intent(out) :: num_flux         ! Output
        real(p2),                   intent(out) :: wsn   

        ! Local Vars
        real(p2), dimension(5) :: qL, qR ! primitive vars reconstructed to face
        real(p2), dimension(5) :: uL, uR ! conservative vars computed from wL and wR
  
        if (accuracy_order == 2) then
            qL = q1 + phi1 * ( gradq1(1,:)*(xm-xc1) + gradq1(2,:)*(ym-yc1) + gradq1(3,:)*(zm-zc1) ) ! gradq <=> gradq (var) 
            qR = q2 + phi2 * ( gradq2(1,:)*(xm-xc2) + gradq2(2,:)*(ym-yc2) + gradq2(3,:)*(zm-zc2) ) ! u <=> q (vars)
        else
            qL = q1 
            qR = q2
        end if

        uL = q2u(qL)
        uR = q2u(qR)

        !------------------------------------------------------------
        !  (1) Roe flux
        !------------------------------------------------------------
        select case(imethod_inv_flux)
        case(IFLUX_ROE)
            call roe(uL,uR,n12, num_flux,wsn)
        ! !------------------------------------------------------------
        ! Other fluxes not yet implemneted.
        !------------------------------------------------------------
        case default

            write(*,*) " Invalid input for inviscid_flux = ", imethod_inv_flux
            write(*,*) " Choose roe or rhll, and try again."
            write(*,*) " ... Stop."
            stop

        end select

    end subroutine interface_flux

end module interface