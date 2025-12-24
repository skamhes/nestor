module interface

    implicit none
    
    contains
    
    subroutine interface_flux(qL, qR, n12, ur2_1, ur2_2, &
                                            num_flux, wsn  )

        use common                 , only : p2

        use config                 , only : method_inv_flux !name of flux pecified as input

        use config                 , only : accuracy_order

        use utils                  , only : imethod_inv_flux, IFLUX_ROE

        use solution               , only : q2u

        use inviscid_flux          , only : roe, roe_lm_w
        implicit none

        ! Inputs
        real(p2), dimension(5),     intent(in) :: qL, qR            ! Primative vars (prim for low mach)
        real(p2), dimension(3),     intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                   intent(in)  :: ur2_1, ur2_2     ! Reference velocity for roe_lm_ws

        ! Output
        real(p2), dimension(5),     intent(out) :: num_flux         ! Output
        real(p2),                   intent(out) :: wsn   

        ! Optional
        

        ! Local Vars
        real(p2), dimension(5) :: uL, uR, num_flux1 ! conservative vars computed from wL and wR


        !------------------------------------------------------------
        !  (1) Roe flux
        !------------------------------------------------------------
        select case(imethod_inv_flux)
        case(IFLUX_ROE)
            call roe(qL,qR,n12,num_flux,wsn)

        elseif(trim(method_inv_flux)=="roe_lm_w") then
            call roe_lm_w(uL,uR,ur2_1,ur2_2,n12, num_flux,wsn)
        !------------------------------------------------------------
        ! Other fluxes not yet implemneted.
        !------------------------------------------------------------
        case default

            write(*,*) " Invalid input for inviscid_flux = ", imethod_inv_flux
            write(*,*) " Choose roe or rhll, and try again."
            write(*,*) " ... Stop."
            stop

        end select

    end subroutine interface_flux

        ! Reconstruct flow up to a cell face
    subroutine reconstruct_flow(xc,xf,phi,q1,gradq, qf)

        use common , only : p2

        implicit none

        real(p2), dimension(3),   intent(in ) :: xc, xf   ! Cell and face centers
        real(p2),                 intent(in ) :: phi      ! Gradient limiter
        real(p2), dimension(5),   intent(in ) :: q1       ! Cell average (1st order) flow values
        real(p2), dimension(3,5), intent(in ) :: gradq    ! Cell gradient
        
        real(p2), dimension(5),   intent(out) :: qf       ! Reconstructed face variables
        
        real(p2), dimension(3)                :: dx

        dx = xf - xc

        qf = q1 + phi * matmul(dx,gradq) 
        !qf = q1 + phi * ( gradq(1,:)*(xf(1)-xc(1)) + gradq(2,:)*(xf(2)-xc(2)) + gradq(3,:)*(xf(3)-xc(3)) )
        ! The second thing is the same just (presumably) slower. (Godbolt says yes!)
    end subroutine reconstruct_flow

end module interface