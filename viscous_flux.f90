module viscous_flux

    implicit none

    public :: visc_flux_internal
    public :: visc_flux_boundary
    public :: compute_visc_num_flux

    private 

    contains

    subroutine visc_flux_internal(q1,q2,gradq1,gradq2,n12,xc1,yc1,zc1,xc2,yc2,zc2, num_flux)

        ! Face gradient terms computed using EQ. 14 in https://doi.org/10.2514/2.689 

        use common                  , only : p2, half, one, zero, three_half, two_third, four_third

        use solution                , only : gammamo, nq, ndim, T_inf ! w2u, nq
        
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
        real(p2)                     :: magds
        real(p2), dimension(nq)      :: qL, qR
        
        integer                      :: ivar

        ! Calculate the face gradients
        ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
        magds = dot_product(ds,ds)
        dsds2 = ds/magds ! ds(:)/ds**2
        magds = sqrt(magds)


        ! write(*,'(a40,3es13.5)') " numerical (uncorrected) face gradx", half * (gradq1(:,2) + gradq2(:,2))
        ! Equation 14
        do ivar = 1,nq
            gradq_face(:,ivar) = half * (gradq1(:,ivar) + gradq2(:,ivar))
            ! gradq_face(:,ivar) = gradq_face(:,ivar) + ( (q2(ivar) - q1(ivar)) - dot_product(gradq_face(:,ivar),ds)) * dsds2
        end do

        ! write(*,'(a40,3es13.5)') " numerical (corrected) face gradx", gradq_face(:,2)
        ! This subroutine only handles computing the interface gradient.
        ! Once we have it we call the internal function
        
        ! We need the reconstructed face values
        ! do ivar = 1,nq
        !     qL(ivar) = q1(ivar) + half * dot_product(gradq1(:,ivar),ds(:))
        !     qR(ivar) = q2(ivar) + half * dot_product(gradq2(:,ivar),ds(:))
        ! end do
        ! ds    = ds / magds
        ! call viscous_alpha_ddt(qL,qR,gradq1,gradq2, n12, ds,magds, num_flux)

        ! write(*,*)
        ! write(*,'(a15,5es13.5)') " alpha_num_flux: ", num_flux
        call compute_visc_num_flux(q1,q2,gradq_face,n12,num_flux)
        ! write(*,'(a15,5es13.5)') " old_num_flux: ", num_flux
        ! write(*,*)
    end subroutine visc_flux_internal

    subroutine visc_flux_boundary(q1,qb,face_gradient,n12,xc1,yc1,zc1,xf2,yf2,zf2,num_flux)

        use common                  , only : p2, half, one, zero, three_half, two_third, four_third

        use solution                , only : gammamo, nq, ndim, T_inf ! w2u, nq
        
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
        real(p2)                     :: magds
        real(p2), dimension(nq)      :: qL, qR
        
        integer                      :: ivar

        gradq_face = face_gradient

        ! Calculate the face gradients
        ds = (/xf2-xc1, yf2-yc1, zf2-zc1/) ! vector pointing from center of cell 1 to cell 2
        
        dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds(:)/ds**2
        
        ! Equation 14
        do ivar = 1,nq
            ! gradq_face(:,ivar) = gradq_face(:,ivar) + ( (qb(ivar) - q1(ivar)) - dot_product(gradq_face(:,ivar),ds)) * dsds2
        end do


        ! This is just a wrapper function since we already have the interface gradient computed.
        call compute_visc_num_flux(q1,qb,gradq_face,n12,num_flux)
        ! ds = ds * 2.0_p2
        ! do ivar = 1,nq
        !     qL(ivar) = q1(ivar) + half * dot_product(face_gradient(:,ivar),ds(:))
        !     qR(ivar) = qb(ivar) + half * dot_product(face_gradient(:,ivar),ds(:))
        ! end do
        ! magds = dot_product(ds,ds)
        ! magds = sqrt(magds)
        ! ds = ds / magds
        ! call viscous_alpha_ddt(qL,qR,face_gradient,face_gradient, n12, ds,magds, num_flux)
    end subroutine visc_flux_boundary

    subroutine compute_visc_num_flux(q1,q2,interface_grad,n12,num_flux)
        use common                  , only : p2, half, one, zero, three_half, two_third, four_third, ix, iy, iz

        use solution                , only : gammamo, nq, ndim, T_inf ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        implicit none 

        real(p2), dimension(nq),      intent(in)    :: q1, q2
        real(p2), dimension(ndim,nq), intent(in)    :: interface_grad
        real(p2), dimension(ndim),    intent(in)    :: n12
        real(p2), dimension(nq),      intent(out)   :: num_flux

        ! Local Vars
        real(p2)                     :: mu
        real(p2)                     :: u, v, w, T
        real(p2)                     :: C0
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
        C0= sutherland_constant/reference_temp
        mu =  M_inf/Re_inf * (one + C0) / (T + C0)*T**(three_half)

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
       
        tauxx =  mu * (four_third*grad_u(ix) - two_third*grad_v(iy) - two_third*grad_w(iz))
        tauyy =  mu * (four_third*grad_v(iy) - two_third*grad_u(ix) - two_third*grad_w(iz))
        tauzz =  mu * (four_third*grad_w(iz) - two_third*grad_u(ix) - two_third*grad_v(iy))

        tauxy =  mu*(grad_u(2) + grad_v(1))
        tauxz =  mu*(grad_u(3) + grad_w(1))
        tauyz =  mu*(grad_v(3) + grad_w(2))
    
        tauyx = tauxy
        tauzx = tauxz
        tauzy = tauyz
    
        ! Heat fluxes: q = - mu*grad(T)/(Prandtl*(gamma-1))
    
        qx = - mu*grad_T(ix)/(pr*(gammamo))
        qy = - mu*grad_T(iy)/(pr*(gammamo))
        qz = - mu*grad_T(iz)/(pr*(gammamo))

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

    !********************************************************************************
!* -- 3D Viscous Flux Function (alpha-damping scheme) and Jacobian ---
!*
!* This subroutine computes an approximate numerical flux for the viscous part of the
!* projected Navier-Stokes flux in the direction njk=[nx,ny,nz].
!*
!* Conservative form of the Navier-Stokes equations:
!*
!*     dU/dt + dF/dx + dG/dy = 0
!*
!* where
!*
!*    F = Fi + Fv   (Fi = inviscid part, Fv = viscous part)
!*    G = Gi + Gv   (Gi = inviscid part, Gv = viscous part)
!*
!* Viscous flux is given by
!*
!*   Fvn = Fv*nx + Gv*ny = |        0      |
!*                         | -  tauxn      |
!*                         | -  tauyn      |
!*                         | -  tauzn      |
!*                         | -  taunv + qn |
!*
!* where 'qn' is the normal heat flux, and
!*
!*      tauxn = tauxx*nx + tauxy*ny + tauxz*nz,
!*      tauyn = tauyx*nx + tauyy*ny + tauyz*nz,
!*      tauzn = tauzx*nx + tauzy*ny + tauzz*nz,
!*      taunv = tauxn*u  + tauyn*v  + tauzn*w.
!*
!* Traditional approach is to individually compute all the quantities in 
!* the viscous flux above at the interface, and then directly evaluate the
!* viscous flux with them. The most critical part is the computation of the
!* interface gradients which are required to compute the viscous stresses
!* and heat fluxes. Gradients must be defined to ensure consistency, design
!* accuracy, robustness, and h-ellipticity.
!*
!* Accurate and robust formulas for the interface gradient are given in
!* Nishikawa-AIAA2010-5093:
!*
!* http://ossanworld.com/hiroakinishikawa/My_papers/nishikawa_AIAA-2010-5093.pdf
!*
!* We employ the one proposed for finite-volume schemes, which has been proved
!* to be robust for unstructured-grid simulations. The value of the damping
!* coefficient, alpha, is taken to be 4/3
!* which corresponds to 4-th order accurate diffusion scheme in 1D.
!* The scheme is a generalized average-least-squares scheme, which includes
!* widely-used schemes (edge-normal and Mathur-Murthy/face-tangent schemes).
!*
!* ------------------------------------------------------------------------------
!*  Input: ucL(1:5) =  Left state (rhoL, rhoL*uL, rhoL*vL, rhoL*wL, rhoL*EL)
!*         ucR(1:5) = Right state (rhoR, rhoR*uR, rhoR*vR, rhoR*wR, rhoR*ER)
!*         njk(1:3) = Unit face vector
!*         ejk(1:3) = Unit edge-vector
!*         mag_ejk  = Magnitude of the edge vector
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | . ejk (edge vector from left to right)
!*              |.
!*       -------x-------- Face
!*             .
!*            .
!*           .
!*          o Left data point
!*
!*  - mag_ejk is the length between the two data points.
!*
!*
!* Output: numerical_flux(1:5) = Numerical viscous flux (alpha-damping scheme),
!*                               and its derivatives as "numerical_flux%df".
!* ------------------------------------------------------------------------------
!*
!* NOTE: This subsroutine computes an approximate version of the alpha-damping
!*       viscous flux (edge-terms-only scheme) which ignores average LSQ gradients.
!*       The numerical flux is therefore not necessarily consistent.
!*       It is consistent only when the skewness parameter = 1.0.
!*       It is OK because this subroutine is designed for computing approximate
!*       derivative of the alpha-damping scheme.
!*
!* Note: You can input the temperature gradient instead of density and pressure
!*       gradients if you have it. I think many practical codes do so.
!*
!* Note: This function is a ddt-version, which means that each variable carries
!*       its derivatives, and that the resulting flux "numerical_flux" will have
!*       its derivatives in "numerical_flux%df". You can convert it to a real-value
!*       version by real(p2) -> real(p2).
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Katate Masatsuka, December 2012. http://www.cfdbooks.com
!********************************************************************************
subroutine viscous_alpha_ddt(qcL,qcR,gradqL,gradqR, njk,ejk,mag_ejk, numerical_flux)

    use common , only : p2, zero, half, two_third, one, four_third, three, ix, iy, iz

    use config , only : sutherland_constant, pr, Re_inf, M_inf

    use solution , only : gamma, T_inf, ip, iu, iv, iw, iT

    implicit none

    !Input
    real(p2), dimension(5)  , intent( in) :: qcL     !Left state (conservative)
    real(p2), dimension(5)  , intent( in) :: qcR     !Right state (conservative)
    real(p2), dimension(3,5), intent( in) :: gradqL  !Left gradient (primitive)
    real(p2), dimension(3,5), intent( in) :: gradqR  !Right gradient (primitive)
    real(p2), dimension(3)  , intent( in) :: njk     !Unit directed area vector
    real(p2), dimension(3)  , intent( in) :: ejk     !Unit edge vector
    real(p2),                 intent( in) :: mag_ejk !Magnitude of the edge vector

    !Output
    real(p2), dimension(5),   intent(out) :: numerical_flux !Numerical viscous flux

    real(p2) ::          C   !Parameter for Sutherland's law
    real(p2) ::    Prandtl    !Prandtl number

    !Local variables
    real(p2) :: alpha !Damping coefficient (see Nishikawa AIAA2010-5093)
    real(p2) :: Lr    !Length scale        (see Nishikawa AIAA2010-5093)

    real(p2) ::   uL, uR            ! x-velocity  (Left and Right states)
    real(p2) ::   vL, vR            ! y-velocity  (Left and Right states)
    real(p2) ::   wL, wR            ! z-velocity  (Left and Right states)
    real(p2) :: rhoL, rhoR          ! Density     (Left and Right states)
    real(p2) ::   presL, presR            ! Pressure    (Left and Right states)
    real(p2) ::   TL, TR            ! Temperature (Left and Right states)

    real(p2) :: tauxx, tauyy, tauzz !Viscous stresses: diagonal compontens
    real(p2) :: tauxy, tauyz, tauzx !Viscous stresses: off-diagonal components
    real(p2) :: tauyx, tauzy, tauxz !Viscous stresses: same as above by symmetry
    real(p2) :: qx, qy, qz          !Heat flux components
    real(p2) :: tauxn, tauyn, tauzn !Normal stresses
    real(p2) :: qn                  !Normal heat flux

    real(p2) :: u, v, w, T, mu                         !Interface quantities
    real(p2), dimension(3) :: grad_u, grad_v, grad_w   !Interface gradients of velocities
    real(p2), dimension(3) :: grad_rho, grad_p, grad_T !Interface gradients of rho, p, and T

    real(p2), dimension(3) :: grad_uL, grad_vL, grad_wL, grad_rL, grad_pL, grad_TL
    real(p2), dimension(3) :: grad_uR, grad_vR, grad_wR, grad_rR, grad_pR, grad_TR

    real(p2) :: rho, a2   !Interface values for density and (speed of sound)^2

    C = sutherland_constant
    Prandtl = pr
    

    ! Left and right states and gradients:

    ! rhoL = ucL(1)
    presL = qcL(ip)
    uL = qcL(iu)
    vL = qcL(iv)
    wL = qcL(iw)
    TL = qcL(iT)
    
    ! grad_rL = gradwL(1,:)
    grad_pL = gradqL(:,ip)
    grad_uL = gradqL(:,iu)
    grad_vL = gradqL(:,iv)
    grad_wL = gradqL(:,iw)
    grad_TL = gradqL(:,iT)

    ! rhoR = ucR(1)
    presR = qcL(ip)
    uR = qcR(iu)
    vR = qcR(iv)
    wR = qcR(iw)
    TR = qcR(iT)

    ! grad_rR = gradwR(1,:)
    grad_pR = gradqR(:,ip)
    grad_uR = gradqR(:,iu)
    grad_vR = gradqR(:,iv)
    grad_wR = gradqR(:,iw)
    grad_TR = gradqR(:,iT)

    ! Arithmetic averages of velocities and temperature.

    u = half*(uL + uR)
    v = half*(vL + vR)
    w = half*(wL + wR)
    T = half*(TL + TR)

    ! Sutherland's law in the nondimensional form.
    ! Note: The factor, M_inf/Re_inf, comes from nondimensionalization.

    mu = (one+C/T_inf)/(T+C/T_inf)*T**(three*half) * M_inf/Re_inf

    ! Damping coefficient, alpha:
    !  (1)alpha=1   gives the central-difference formula in 1D.
    !  (2)alpha=4/3 corresponds to the 4th-order diffusion scheme in 1D.

    alpha = one
    ! alpha = four_third

    ! Lr = Length scale involving the skewness measure (see Nishikawa AIAA2010-5093)
    ! This is the key quantity for robust and accurate computations on skewed grids.

    Lr = half*abs( njk(ix)*ejk(ix) + njk(iy)*ejk(iy) + njk(iz)*ejk(iz) ) * mag_ejk

    ! Interface gradients from the derived diffusion scheme (Nishikawa-AIAA2010-5093).
    ! The second term is the damping term.

    grad_u = half*( (grad_uR + grad_uL) + alpha/Lr*(uR-uL)*njk )
    grad_v = half*( (grad_vR + grad_vL) + alpha/Lr*(vR-vL)*njk )
    grad_w = half*( (grad_wR + grad_wL) + alpha/Lr*(wR-wL)*njk )
    grad_T = half*( (grad_TR + grad_TL) + alpha/Lr*(TR-TL)*njk )

    ! The temperature gradient is computed from the interface density and the pressure
    ! gradients. 
    ! Note: T = gamma*p/rho -> grad(T) = gamma*grad(p)/rho - (gamma*p/rho^2)*grad(rho)

    ! rho = half*(rhoR + rhoL)
    ! a2 = gamma*half*(pR + pL)/rho

    ! grad_rho = half*( (grad_rR + grad_rL) + alpha/Lr*(rhoR-rhoL)*njk )
    ! grad_p   = half*( (grad_pR + grad_pL) + alpha/Lr*(  pR-pL  )*njk )

    ! grad_T   = ( gamma*grad_p - a2*grad_rho) /rho

    ! Interface gradients have been computed: grad(u), grad(v), grad(w), grad(T).
    ! We now evaluate the physical viscous flux with them.

    ! Viscous stresses (Stokes' hypothesis is assumed)

    tauxx =  mu*(four_third*grad_u(ix) - two_third*grad_v(iy) - two_third*grad_w(iz))
    tauyy =  mu*(four_third*grad_v(iy) - two_third*grad_u(ix) - two_third*grad_w(iz))
    tauzz =  mu*(four_third*grad_w(iz) - two_third*grad_u(ix) - two_third*grad_v(iy))

    tauxy =  mu*(grad_u(iy) + grad_v(ix))
    tauxz =  mu*(grad_u(iz) + grad_w(ix))
    tauyz =  mu*(grad_v(iz) + grad_w(iy))

    tauyx = tauxy
    tauzx = tauxz
    tauzy = tauyz

    ! Heat fluxes: q = - mu*grad(T)/(Prandtl*(gamma-1))

    qx = - mu*grad_T(ix)/(Prandtl*(gamma-one))
    qy = - mu*grad_T(iy)/(Prandtl*(gamma-one))
    qz = - mu*grad_T(iz)/(Prandtl*(gamma-one))

    ! Normal components

    tauxn = tauxx*njk(ix) + tauxy*njk(iy) + tauxz*njk(iz)
    tauyn = tauyx*njk(ix) + tauyy*njk(iy) + tauyz*njk(iz)
    tauzn = tauzx*njk(ix) + tauzy*njk(iy) + tauzz*njk(iz)
    qn    = qx*njk(ix)    + qy*njk(iy)    + qz*njk(iz)

    ! Evaluate the viscous flux at the interface

    numerical_flux(1) =   zero
    numerical_flux(2) = - tauxn
    numerical_flux(3) = - tauyn
    numerical_flux(4) = - tauzn
    numerical_flux(5) = - (tauxn*u + tauyn*v + tauzn*w) + qn

    ! Normal max wave speed
    !    wsn = alpha*(mu/rho*gamma/Prandtl)/Lr

end subroutine viscous_alpha_ddt
!--------------------------------------------------------------------------------

end module viscous_flux