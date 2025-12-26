module viscous_flux

    implicit none

    public :: visc_flux_internal
    public :: visc_flux_boundary

    private 

    contains

    subroutine visc_flux_internal(q1,q2,gradq1,gradq2,trb1,trb2,n12,xc1,yc1,zc1,xc2,yc2,zc2, num_flux)

        ! Face gradient terms computed using EQ. 14 in https://doi.org/10.2514/2.689 

        use common                  , only : p2, half

        use solution_vars           , only : nq, ndim ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use turb                    , only : nturb

        implicit none

        real(p2), dimension(nq),      intent(in) :: q1, q2
        real(p2), dimension(ndim,nq), intent(in) :: gradq1, gradq2
        real(p2), dimension(nturb),   intent(in) :: trb1, trb2
        real(p2), dimension(ndim),    intent(in) :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                     intent(in) :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                     intent(in) :: xc2, yc2, zc2     ! Right cell centroid
        real(p2), dimension(nq),      INTENT(OUT):: num_flux

        ! Local Vars
        real(p2), dimension(ndim,nq) :: gradq_face
        real(p2), dimension(ndim)    :: ds,  dsds2
        real(p2)                     :: magds
        real(p2), dimension(nq)      :: correction
        
        integer                      :: ivar

        ! Calculate the face gradients
         ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
        magds = dot_product(ds,ds)
        dsds2 = ds/magds ! ds(:)/ds**2
        ! magds = sqrt(magds)

        ! Equation 14 (optimized to save ~1e-06 second per call) But like... this gets called a lot man.
         gradq_face = half*(gradq1+gradq2)
         correction = matmul(ds,gradq_face) - (q2-q1)
         do ivar = 2,5
             gradq_face(:,ivar) = gradq_face(:,ivar) - correction(ivar) * dsds2
         end do

        ! This subroutine only handles computing the interface gradient.
        ! Once we have it we call the internal function
        call compute_visc_num_flux(q1,q2,trb1,trb2,gradq_face,n12,num_flux)

    end subroutine visc_flux_internal

    subroutine visc_flux_boundary(q1,qb,trb1,trb2,face_gradient,n12,xc1,yc1,zc1,xc2,yc2,zc2,num_flux)

        use common                  , only : p2, half

        use solution_vars           , only : nq, ndim ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use turb                    , only : nturb

        implicit none

        real(p2), dimension(nq),      intent(in)    :: q1, qb
        real(p2), dimension(nturb),   intent(in)    :: trb1, trb2
        real(p2), dimension(ndim,nq), intent(in)    :: face_gradient     ! Grad at bound interface computed using avg face's vgrad
        real(p2), dimension(ndim),    intent(in)    :: n12               ! Normalized face vector
        real(p2),                     intent(in)    :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                     intent(in)    :: xc2, yc2, zc2     ! Boundary face centroid
        real(p2), dimension(nq),      intent(out)   :: num_flux


        ! Local Vars
        real(p2), dimension(ndim,nq) :: gradq_face
        real(p2), dimension(ndim)    :: ds,  dsds2
        real(p2)                     :: magds
        real(p2), dimension(nq)      :: correction

        integer                      :: ivar

        ! gradq_face = face_gradient

        ! Calculate the face gradients
        ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
        magds = dot_product(ds,ds)
        dsds2 = ds/magds ! ds(:)/ds**2
        ! magds = sqrt(magds)

        correction = matmul(ds,face_gradient) - (qb-q1)
        do ivar = 2,5
            gradq_face(:,ivar) = face_gradient(:,ivar) - correction(ivar) * dsds2
        end do
        ! This is just a wrapper function since we already have the interface gradient computed.
        call compute_visc_num_flux(q1,qb,trb1,trb2,gradq_face,n12,num_flux)

        
    end subroutine visc_flux_boundary

    subroutine compute_visc_num_flux(q1,q2,trb1,trb2,interface_grad,n12,num_flux)
        use common                  , only : p2, half, zero, two_third, four_third, ix, iy, iz

        use solution_vars           , only : gammamo, nq, ndim, iu, iv, iw, iT, C0 ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp, pr_t

        use viscosity               , only : compute_viscosity

        use turb                    , only : nturb, calcmut

        implicit none 

        real(p2), dimension(nq),      intent(in)    :: q1, q2
        real(p2), dimension(nturb),   intent(in)    :: trb1, trb2
        real(p2), dimension(ndim,nq), intent(in)    :: interface_grad
        real(p2), dimension(ndim),    intent(in)    :: n12
        real(p2), dimension(nq),      intent(out)   :: num_flux

        ! Local Vars
        real(p2)                     :: mu_effective, mu, mut
        real(p2), dimension(nturb)   :: trb
        real(p2), dimension(nq)      :: qf
        real(p2)                     :: p, u, v, w, T
        real(p2)                     :: tauxx, tauyy, tauzz !Viscous stresses: diagonal compontens
        real(p2)                     :: tauxy, tauyz, tauzx !Viscous stresses: off-diagonal components
        real(p2)                     :: tauyx, tauzy, tauxz !Viscous stresses: same as above by symmetry
        real(p2)                     :: qx, qy, qz          !Heat flux components
        real(p2)                     :: tauxn, tauyn, tauzn !Normal stresses
        real(p2)                     :: qn                  !Normal heat flux
        real(p2), dimension(3)       :: grad_u, grad_v, grad_w   !Interface gradients of velocities
        real(p2), dimension(3)       :: grad_T
        
        qf = half * (q1 + q2)

        p = qf(1)  ! p at the face
        u = qf(2)  ! u at the face
        v = qf(3)  ! v at the face
        w = qf(4)  ! w at the face
        T = qf(5)  ! T at the face
        
        mu = compute_viscosity(T)

        trb = half * (trb1 + trb2)

        mut = calcmut(qf,mu,trb)

        mu_effective = mu + mut

#ifdef NANCHECK
        if (isnan(mu_effective)) then 
            write (*,*) "nan value present - press [Enter] to continue"
            read(unit=*,fmt=*)
        end if
#endif
        ! Interface values
        grad_u = interface_grad(:,iu)
        grad_v = interface_grad(:,iv)
        grad_w = interface_grad(:,iw)
        grad_T = interface_grad(:,iT)

        ! Viscous stresses (Stokes' hypothesis is assumed)
       
        tauxx =  mu_effective*(four_third*grad_u(ix) - two_third*grad_v(iy) - two_third*grad_w(iz))
        tauyy =  mu_effective*(four_third*grad_v(iy) - two_third*grad_u(ix) - two_third*grad_w(iz))
        tauzz =  mu_effective*(four_third*grad_w(iz) - two_third*grad_u(ix) - two_third*grad_v(iy))
    
        tauxy =  mu_effective*(grad_u(iy) + grad_v(ix))
        tauxz =  mu_effective*(grad_u(iz) + grad_w(ix))
        tauyz =  mu_effective*(grad_v(iz) + grad_w(iy))
    
        tauyx = tauxy
        tauzx = tauxz
        tauzy = tauyz
    
        ! Heat fluxes: q = - mu*grad(T)/(Prandtl*(gamma-1))
    
        qx = - ( mu/(pr*gammamo) + mut/(pr_t*gammamo) ) * grad_T(1)
        qy = - ( mu/(pr*gammamo) + mut/(pr_t*gammamo) ) * grad_T(2)
        qz = - ( mu/(pr*gammamo) + mut/(pr_t*gammamo) ) * grad_T(3)

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