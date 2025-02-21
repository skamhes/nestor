module viscous_flux

    implicit none

    public :: visc_flux_internal
    public :: visc_flux_boundary

    private 

    contains

    subroutine visc_flux_internal(q1,q2,gradq1,gradq2,n12,xc1,yc1,zc1,xc2,yc2,zc2, num_flux)

        ! Face gradient terms computed using EQ. 14 in https://doi.org/10.2514/2.689 

        use common                  , only : p2, half

        use solution_vars           , only : nq, ndim ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        implicit none

        real(p2), dimension(nq),      intent(in) :: q1, q2
        real(p2), dimension(ndim,nq), intent(in) :: gradq1, gradq2
        real(p2), dimension(ndim),    intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                     intent(in) :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                     intent(in) :: xc2, yc2, zc2     ! Right cell centroid
        real(p2), dimension(nq),      INTENT(OUT):: num_flux

        ! Local Vars
        real(p2), dimension(ndim,nq) :: gradq_face
        real(p2), dimension(ndim)    :: ds,  dsds2
        
        integer                      :: ivar

        ! Calculate the face gradients
        ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
        dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds(:)/ds**2

        ! Equation 14
        do ivar = 1,nq
            gradq_face(:,ivar) = half * (gradq1(:,ivar) + gradq2(:,ivar))
            gradq_face(:,ivar) = gradq_face(:,ivar) + ( (q2(ivar) - q1(ivar)) - dot_product(gradq_face(:,ivar),ds)) * dsds2
        end do

        ! This subroutine only handles computing the interface gradient.
        ! Once we have it we call the internal function
        call compute_visc_num_flux(q1,q2,gradq_face,n12,num_flux)

    end subroutine visc_flux_internal

    subroutine visc_flux_boundary(q1,qb,face_gradient,n12,xc1,yc1,zc1,xf2,yf2,zf2,num_flux)

        use common                  , only : p2, half

        use solution_vars           , only : nq, ndim ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        implicit none

        real(p2), dimension(nq),      intent(in)    :: q1, qb
        real(p2), dimension(ndim,nq), intent(in)    :: face_gradient     ! Grad at bound interface computed using avg face's vgrad
        real(p2), dimension(ndim),    intent(in)    :: n12               ! Normalized face vector
        real(p2),                     intent(in)    :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                     intent(in)    :: xf2, yf2, zf2     ! Boundary face centroid
        real(p2), dimension(nq),      intent(out)   :: num_flux


        ! Local Vars
        real(p2), dimension(ndim,nq) :: gradq_face
        real(p2), dimension(ndim)    :: ds,  dsds2
        
        integer                      :: ivar

        gradq_face = face_gradient

        ! Calculate the face gradients
        ds = (/xf2-xc1, yf2-yc1, zf2-zc1/) ! vector pointing from center of cell 1 to cell 2
        dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds(:)/ds**2

        ! Equation 14
        do ivar = 1,nq
            gradq_face(:,ivar) = gradq_face(:,ivar) + ( half * (qb(ivar) - q1(ivar)) - dot_product(gradq_face(:,ivar),ds)) * dsds2
        end do


        ! This is just a wrapper function since we already have the interface gradient computed.
        call compute_visc_num_flux(q1,qb,gradq_face,n12,num_flux)

        
    end subroutine visc_flux_boundary

    subroutine compute_visc_num_flux(q1,q2,interface_grad,n12,num_flux)
        use common                  , only : p2, half, zero, two_third, four_third

        use solution_vars           , only : gammamo, nq, ndim ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use viscosity               , only : compute_viscosity
        implicit none 

        real(p2), dimension(nq),      intent(in)    :: q1, q2
        real(p2), dimension(ndim,nq), intent(in)    :: interface_grad
        real(p2), dimension(ndim),    intent(in)    :: n12
        real(p2), dimension(nq),      intent(out)   :: num_flux

        ! Local Vars
        real(p2)                     :: mu
        real(p2)                     :: u, v, w, T
        real(p2)                     :: tauxx, tauyy, tauzz !Viscous stresses: diagonal compontens
        real(p2)                     :: tauxy, tauyz, tauzx !Viscous stresses: off-diagonal components
        real(p2)                     :: tauyx, tauzy, tauxz !Viscous stresses: same as above by symmetry
        real(p2)                     :: qx, qy, qz          !Heat flux components
        real(p2)                     :: tauxn, tauyn, tauzn !Normal stresses
        real(p2)                     :: qn                  !Normal heat flux
        real(p2), dimension(3)       :: grad_u, grad_v, grad_w   !Interface gradients of velocities
        real(p2), dimension(3)       :: grad_T
        
        u = half * (q1(2)  + q2(2) ) ! u at the face
        v = half * (q1(3)  + q2(3) ) ! v at the face
        w = half * (q1(4)  + q2(4) ) ! w at the face
        T = half * (q1(nq) + q2(nq)) ! T at the face
        
        mu = compute_viscosity(T)

        ! get_viscosity = scaling_factor * ( (one + ( C_0/Freestream_Temp ) )/(T + ( C_0/Freestream_Temp )) ) ** 1.5_p2
        if (isnan(mu)) then 
            write (*,*) "nan value present - press [Enter] to continue"
            read(unit=*,fmt=*)
        end if

        ! Interface values
        grad_u = interface_grad(:,2)
        grad_v = interface_grad(:,3)
        grad_w = interface_grad(:,4)
        grad_T = interface_grad(:,5)

        ! Viscous stresses (Stokes' hypothesis is assumed)
       
        tauxx =  mu*(four_third*grad_u(1) - two_third*grad_v(2) - two_third*grad_w(3))
        tauyy =  mu*(four_third*grad_v(2) - two_third*grad_u(1) - two_third*grad_w(3))
        tauzz =  mu*(four_third*grad_w(3) - two_third*grad_u(1) - two_third*grad_v(2))
    
        tauxy =  mu*(grad_u(2) + grad_v(1))
        tauxz =  mu*(grad_u(3) + grad_w(1))
        tauyz =  mu*(grad_v(3) + grad_w(2))
    
        tauyx = tauxy
        tauzx = tauxz
        tauzy = tauyz
    
        ! Heat fluxes: q = - mu*grad(T)/(Prandtl*(gamma-1))
    
        qx = - mu*grad_T(1)/(pr*(gammamo))
        qy = - mu*grad_T(2)/(pr*(gammamo))
        qz = - mu*grad_T(3)/(pr*(gammamo))

        tauxn = tauxx*n12(1) + tauxy*n12(2) + tauxz*n12(3)
        tauyn = tauyx*n12(1) + tauyy*n12(2) + tauyz*n12(3)
        tauzn = tauzx*n12(1) + tauzy*n12(2) + tauzz*n12(3)
        qn    =    qx*n12(1) +    qy*n12(2) +    qz*n12(3)

        num_flux(1) =   zero
        num_flux(2) = - tauxn
        num_flux(3) = - tauyn
        num_flux(4) = - tauzn
        num_flux(5) = - (tauxn*u + tauyn*v + tauzn*w) + qn

    end subroutine compute_visc_num_flux

end module viscous_flux