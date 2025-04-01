module inviscid_flux

    implicit none

    public

    contains

!********************************************************************************
    !* -- 3D Roe's Flux Function and Jacobian --
    !*
    !* NOTE: This version does not use any tangent vector.
    !*       See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
    !*
    !* This subroutine computes the Roe flux for the Euler equations
    !* in the direction, njk=[nx,ny,nz].
    !*
    !* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
    !* Schemes, Journal of Computational Physics, 43, pp. 357-372.
    !*
    !* Conservative form of the Euler equations:
    !*
    !*     dU/dt + dF/dx + dG/dy + dH/dz = 0
    !*
    !* This subroutine computes the numerical flux for the flux in the direction,
    !* njk=[nx,ny,nz]:
    !*
    !*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
    !*                               | rho*qn*u + p*nx |
    !*                               | rho*qn*v + p*ny |
    !*                               | rho*qn*w + p*nz |
    !*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
    !*
    !* The Roe flux is implemented in the following form:
    !*
    !*   Numerical flux = 1/2 [ Fn(UR) + Fn(UL) - |An|dU ], 
    !*
    !*  where
    !*
    !*    An = dFn/dU,  |An| = R|Lambda|L, dU = UR - UL.
    !*
    !* The dissipation term, |An|dU, is actually computed as
    !*
    !*     sum_{k=1,4} |lambda_k| * (LdU)_k * r_k,
    !*
    !* where lambda_k is the k-th eigenvalue, (LdU)_k is the k-th wave strength,
    !* and r_k is the k-th right-eigenvector evaluated at the Roe-average state.
    !*
    !* Note: The 4th component is a combined contribution from two shear waves.
    !*       They are combined to eliminate the tangent vectors.
    !*       So, (LdU)_4 is not really a wave strength, and
    !*       r_4 is not really an eigenvector.
    !*       See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
    !*
    !* Note: In the code, the vector of conserative variables are denoted by uc.
    !*
    !* ------------------------------------------------------------------------------
    !*  Input: ucL(1:5) =  Left state (rhoL, rhoL*uL, rhoL*vL, rhoL*wR, rhoL*EL)
    !*         ucR(1:5) = Right state (rhoR, rhoL*uR, rhoL*vR, rhoL*wR, rhoL*ER)
    !*         njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
    !*
    !*           njk
    !*  Face normal ^   o Right data point
    !*              |  .
    !*              | .
    !*              |. 
    !*       -------x-------- Face
    !*             .                 Left and right states are
    !*            .                   1. Values at data points for 1st-order accuracy
    !*           .                    2. Extrapolated values at the face midpoint 'x'
    !*          o Left data point        for 2nd/higher-order accuracy.
    !*
    !*
    !* Output:  num_flux(1:5) = the numerical flux vector
    !*                    wsn = maximum wave speed (eigenvalue)
    !*
    !* ------------------------------------------------------------------------------
    !*
    !* Note: This subroutine has been prepared for an educational purpose.
    !*       It is not at all efficient. Think about how you can optimize it.
    !*       One way to make it efficient is to reduce the number of local variables,
    !*       by re-using temporary variables as many times as possible.
    !*
    !* Note: Please let me know if you find bugs. I'll greatly appreciate it and
    !*       fix the bugs.
    !*
    !* Katate Masatsuka, November 2012. http://www.cfdbooks.com
    !********************************************************************************
    ! 
    !
    !
    ! Note: This is currently unaltered from the edu_euler roe function (other than some formatting changes).  
    ! I will at some point optimie this.  But right now I'm just looking to get a working code...
    subroutine roe(ucL, ucR, njk, num_flux,wsn)

        use solution_vars    , only : gamma
        use config      , only : eig_limiting_factor
       
        implicit none
       
        integer , parameter :: p2 = selected_real_kind(15) ! Double precision
       
        !Input
        real(p2), dimension(5), intent( in) :: ucL !Left  state in conservative variables.
        real(p2), dimension(5), intent( in) :: ucR !Right state in conservative variables.
        real(p2), dimension(3), intent( in) :: njk
       
        !Output
        real(p2), dimension(5)  , intent(out) :: num_flux !Numerical viscous flux
        real(p2),                 intent(out) :: wsn      !Max wave speed
       
        !Some constants
        real(p2) ::  zero = 0.0_p2
        real(p2) ::   one = 1.0_p2
        real(p2) ::   two = 2.0_p2
        real(p2) ::  half = 0.5_p2
       
        !Local variables
        !            L = Left
        !            R = Right
        ! No subscript = Roe average
       
        real(p2) :: nx, ny, nz             ! Normal vector components
        real(p2) :: uL, uR, vL, vR, wL, wR ! Velocity components.
        real(p2) :: rhoL, rhoR, pL, pR     ! Primitive variables.
        real(p2) :: qnL, qnR               ! Normal velocities
        real(p2) :: aL, aR, HL, HR         ! Speed of sound, Total enthalpy
        real(p2), dimension(5)   :: fL     ! Physical flux evaluated at ucL
        real(p2), dimension(5)   :: fR     ! Physical flux evaluated at ucR
       
        real(p2) :: RT                     ! RT = sqrt(rhoR/rhoL)
        real(p2) :: rho,u,v,w,H,a,qn       ! Roe-averages
       
        real(p2) :: drho, dqn, dp          ! Differences in rho, qn, p, e.g., dp=pR-pL
        real(p2), dimension(4) :: LdU      ! Wave strengths = L*(UR-UL)
        real(p2) :: du, dv, dw             ! Velocity differences
        real(p2), dimension(4)   :: ws     ! Wave speeds
        real(p2), dimension(4)   :: dws    ! Width of a parabolic fit for entropy fix
        real(p2), dimension(5,4) :: R      ! Right-eigenvector matrix
        real(p2), dimension(5)   :: diss   ! Dissipation term
       
        integer  :: i
       
        ! Face normal vector (unit vector)
       
        nx = njk(1)
        ny = njk(2)
        nz = njk(3)
       
        !Primitive and other variables.
        
        !  Left state
       
           rhoL = ucL(1)
             uL = ucL(2)/ucL(1)
             vL = ucL(3)/ucL(1)
             wL = ucL(4)/ucL(1)
            qnL = uL*nx + vL*ny + wL*nz
             pL = (gamma-one)*( ucL(5) - half*rhoL*(uL*uL+vL*vL+wL*wL) )
             aL = sqrt(gamma*pL/rhoL)
             HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
       
        !  Right state
       
           rhoR = ucR(1)
             uR = ucR(2)/ucR(1)
             vR = ucR(3)/ucR(1)
             wR = ucR(4)/ucR(1)
            qnR = uR*nx + vR*ny + wR*nz
             pR = (gamma-one)*( ucR(5) - half*rhoR*(uR*uR+vR*vR+wR*wR) )
             aR = sqrt(gamma*pR/rhoR)
             HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)
       
        !Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)
       
        fL(1) = rhoL*qnL
        fL(2) = rhoL*qnL * uL + pL*nx
        fL(3) = rhoL*qnL * vL + pL*ny
        fL(4) = rhoL*qnL * wL + pL*nz
        fL(5) = rhoL*qnL * HL
    
        fR(1) = rhoR*qnR
        fR(2) = rhoR*qnR * uR + pR*nx
        fR(3) = rhoR*qnR * vR + pR*ny
        fR(4) = rhoR*qnR * wR + pR*nz
        fR(5) = rhoR*qnR * HR
       
        !First compute the Roe-averaged quantities
        
        !  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
        !        the Roe-averaged density.
        
           RT = sqrt(rhoR/rhoL)
          rho = RT*rhoL                                        !Roe-averaged density
            u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
            v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
            w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
            H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
            a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
           qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity
       
        !Wave Strengths
       
          drho = rhoR - rhoL !Density difference
            dp =   pR - pL   !Pressure difference
           dqn =  qnR - qnL  !Normal velocity difference
       
        LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
        LdU(2) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
        LdU(3) =  drho - dp/(a*a)            !Entropy wave strength
        LdU(4) = rho                         !Shear wave strength (not really, just a factor)
       
        !Absolute values of the wave Speeds
    
        ws(1) = abs(qn-a) !Left-moving acoustic wave
        ws(2) = abs(qn+a) !Right-moving acoustic wave
        ws(3) = abs(qn)   !Entropy wave
        ws(4) = abs(qn)   !Shear waves
       
        ! Harten's Entropy Fix JCP(1983), 49, pp357-393. This is typically applied
        ! only for the nonlinear fields (k=1 and 3), but here it is applied to all
        ! for robustness, avoiding vanishing wave speeds by making a parabolic fit
        ! near ws = 0 for all waves.
        ! 02-27-2018: The limiting can be too much for the shear wave and entropy wave.
        !             Flat plate calculation shows that applying it to all contaminates
        !             the solution significantly. So, apply only to the nonlinear waves,
        !             or apply very small limiting to entropy and shear waves.
        !
        ! Note: ws(1) and ws(2) are the nonlinear waves.
       
        do i = 1, 4
            dws(i) = eig_limiting_factor(i)*a
            if ( ws(i) < dws(i) ) ws(i) = half * ( ws(i)*ws(i)/dws(i)+dws(i) )
        end do
       
        !Right Eigenvectors
        !Note: Two shear wave components are combined into one, so that tangent vectors
        !      are not required. And that's why there are only 4 vectors here.
        !      See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
        
        ! Left-moving acoustic wave
        R(1,1) = one    
        R(2,1) = u - a*nx
        R(3,1) = v - a*ny
        R(4,1) = w - a*nz
        R(5,1) = H - a*qn
       
        ! Right-moving acoustic wave
        R(1,2) = one
        R(2,2) = u + a*nx
        R(3,2) = v + a*ny
        R(4,2) = w + a*nz
        R(5,2) = H + a*qn
       
        ! Entropy wave
        R(1,3) = one
        R(2,3) = u
        R(3,3) = v 
        R(4,3) = w
        R(5,3) = half*(u*u + v*v + w*w)
       
        ! Two shear wave components combined into one (wave strength incorporated).
        du = uR - uL
        dv = vR - vL
        dw = wR - wL
        R(1,4) = zero
        R(2,4) = du - dqn*nx
        R(3,4) = dv - dqn*ny
        R(4,4) = dw - dqn*nz
        R(5,4) = u*du + v*dv + w*dw - qn*dqn
       
        !Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]
       
        diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
                + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)
       
        ! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]
        ! write(*,*) "Normal:", diss
        num_flux = half * (fL + fR - diss)
       
        ! Max wave speed normal to the face:
                    wsn = abs(qn) + a
       
        if (any(isnan(num_flux(:)))) then 
          write (*,*) "nan value present - press [Enter] to continue"
          read(unit=*,fmt=*)
        end if
    end subroutine roe
    !--------------------------------------------------------------------------------

end module inviscid_flux