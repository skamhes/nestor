module mms

    use common , only : p2, pi

    implicit none

    private
    public fMMS

    ! coefficients
    real(p2), parameter :: arx = 1.5_p2
    real(p2), parameter :: d = 0.25
    real(p2), parameter :: pi2 = 2.0_p2 * pi

    contains

    subroutine fMMS(x,y,z, Q,S, gradW)

        use config , only : M_inf, Pr, sutherland_constant, Re_inf, M_inf, reference_temp

        use utils , only : iturb_type, TURB_INVISCID

        use common , only : half, one

        use solution , only : w2q, gamma, gammamo

        use mms_funcs

        implicit none

        real(p2),                 intent(in)  :: x,y,z ! coordinates
        real(p2), dimension(5),   intent(out) :: Q, S  ! Solution and source terms
        real(p2), dimension(3,5), optional, intent(inout) :: gradW

        ! Coefficients
        real(p2) :: cr0, crs, crx, cry
        real(p2) :: cu0, cus, cux, cuy
        real(p2) :: cv0, cvs, cvx, cvy
        real(p2) :: cp0, cps, cpx, cpy
        real(p2) ::  a,  ax,  ay
        real(p2) ::  b,  bx,  by
        real(p2) ::  k,  kx,  ky
       
        real(p2) :: u2   , u2x   , u2y
        real(p2) :: v2   , v2x   , v2y
        real(p2) :: rho , u , v , p
        real(p2) :: rhox, ux, vx, px
        real(p2) :: rhoy, uy, vy, py

        real(p2) :: au   , aux   , auy
        real(p2) :: av   , avx   , avy
        real(p2) :: bu   , bux   , buy
        real(p2) :: bv   , bvx   , bvy
        real(p2) :: rhoH , rhoHx , rhoHy
        real(p2) :: rhouH, rhouHx, rhouHy
        real(p2) :: rhovH, rhovHx, rhovHy

        real(p2) :: F2, F3, F5
        ! real(p2) :: F2hand, F3hand, F5hand

        ! ! Stress tensor
        ! real(p2) :: txx, txy
        ! real(p2) :: tyx, tyy
        ! ! d(tau) (we only need some terms)
        ! real(p2) :: txx_x, txy_x, tyx_x
        ! real(p2) :: tyy_y, tyx_y, txy_y
        ! real(p2) :: tvx_x, tvy_y

        ! real(p2) :: t, dtdx, dtdy, dtdz ! temp
        ! real(p2) :: dtdxx, dtdxy, dtdyx, dtdyy
        ! real(p2) :: pxr, prx, pyr, pry
        ! real(p2) :: pxrx, pxry, pyrx, pyry
        ! real(p2) :: pxxr, prxx, pyyr, pryy
        ! real(p2) :: pxyr, prxy, pyxr, pryx
        ! real(p2) :: qxx, qyy, qzz

        ! real(p2) :: dummy
        ! real(p2) :: f, fx, fy, g, gx, gy

        ! ! 2nd derivatives
        ! ! u
        ! real(p2) :: uxx, uxy! uxy = uyx
        ! real(p2) :: uyx, uyy
        ! ! v
        ! real(p2) :: vxx, vxy! vxy = vyx
        ! real(p2) :: vyx, vyy
        ! ! p
        ! real(p2) :: pxx, pxy! pxy = pyx
        ! real(p2) :: pyx, pyy
        ! ! rho
        ! real(p2) :: rxx, rxy! rxy = ryx
        ! real(p2) :: ryx, ryy
        

        ! real(p2) :: mu, mudt, mux, muy, muz
        real(p2) :: C0, xmr

        real(p2), dimension(5) :: wtmp
        
       !-----------------------------------------------------------
        ! Constants for the exact solution: c0 + cs*sin(cx*x+cy*y).
        !
        ! Note: Make sure the density and pressure are positive.
        ! Note: These values are passed to the subroutine:
        !         manufactured_sol(c0,cs,cx,cy, nx,ny,x,y),
        !       whcih returns the solution value or derivatives.

        !-----------------------------------------
        ! Density    = cr0 + crs*sin(crx*x+cry*y)

        cr0 =  1.12_p2
        crs =  0.15_p2
        crx =  3.12_p2*pi
        cry =  2.92_p2*pi

        !-----------------------------------------
        ! X-velocity = cu0 + cus*sin(cux*x+cuy*y)

        cu0 =  1.32_p2
        cus =  0.06_p2
        cux =  2.09_p2*pi
        cuy =  3.12_p2*pi

        !-----------------------------------------
        ! Y-velocity = cv0 + cvs*sin(cvx*x+cvy*y)

        cv0 =  1.18_p2
        cvs =  0.03_p2
        cvx =  2.15_p2*pi
        cvy =  3.32_p2*pi

        !-----------------------------------------
        ! Pressure   = cp0 + cps*sin(cpx*x+cpy*y)

        cp0 =  1.62_p2
        cps =  0.31_p2
        cpx =  3.79_p2*pi
        cpy =  2.98_p2*pi

        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! Part I: Compute w = [rho,u,v,p] and grad(w).
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------

        !------------------------------------------------------------------------
        ! rho: Density and its 1st derivatives

        rho = manufactured_sol(cr0,crs,crx,cry, 0,0,x,y)
        rhox = manufactured_sol(cr0,crs,crx,cry, 1,0,x,y)
        rhoy = manufactured_sol(cr0,crs,crx,cry, 0,1,x,y)
 
        !------------------------------------------------------------------------
        ! u: x-velocity and its 1st derivatives
 
        u  = manufactured_sol(cu0,cus,cux,cuy, 0,0,x,y)
        ux = manufactured_sol(cu0,cus,cux,cuy, 1,0,x,y)
        uy = manufactured_sol(cu0,cus,cux,cuy, 0,1,x,y)
 
        !------------------------------------------------------------------------
        ! v: y-velocity and its 1st derivatives
 
        v  = manufactured_sol(cv0,cvs,cvx,cvy, 0,0,x,y)
        vx = manufactured_sol(cv0,cvs,cvx,cvy, 1,0,x,y)
        vy = manufactured_sol(cv0,cvs,cvx,cvy, 0,1,x,y)
 
        !------------------------------------------------------------------------
        ! p: pressure and its 1st derivatives
 
        p  = manufactured_sol(cp0,cps,cpx,cpy, 0,0,x,y)
        px = manufactured_sol(cp0,cps,cpx,cpy, 1,0,x,y)
        py = manufactured_sol(cp0,cps,cpx,cpy, 0,1,x,y)

        wtmp = (/ rho, u, v, 0.0_p2, p/)
        Q = w2q(wtmp)

        if (present(gradW)) then
            gradW = 0.0_p2
            gradW(:,1) = (/rhox, rhoy, 0.0_p2/)
            gradW(:,2) = (/ux, uy, 0.0_p2/)
            gradW(:,3) = (/vx, vy, 0.0_p2/)
            gradW(:,4) = 0.0_p2
            gradW(:,5) = (/px, py, 0.0_p2/)
        endif

        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! Part II: Compute the forcing terms.
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------

        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        !  Inviscid terms
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------

        ! The subroutine 'derivatives_ab" computes the product and derivatives of
        ! two variables (a,b): a*b, ax*b+a*bx, ay*b+a*by.


        ! Derivatives of u^2
        call derivatives_ab(u,ux,uy, u,ux,uy,  u2,u2x,u2y)

        ! Derivatives of v^2
        call derivatives_ab(v,vx,vy, v,vx,vy,  v2,v2x,v2y)
        
        ! Derivatives of k=(u^2+v^2)/2
        
        k     = half*(u*u  + v*v)
        kx    = half*(u2x   + v2x)
        ky    = half*(u2y   + v2y)
    
        ! Derivatives of rho*k = rho*(u^2+v^2)/2
        call derivatives_ab(rho,rhox,rhoy, k,kx,ky,  a,ax,ay) !a=rho*(u^2+v^2)/2
        
        ! Derivatives of rho*H = gamma/(gamma-1)*p + rho*k
        rhoH    = gamma/(gamma-one)*p    + a
        rhoHx   = gamma/(gamma-one)*px   + ax
        rhoHy   = gamma/(gamma-one)*py   + ay
        
        !-----------------------------------------------------------------------------
        
        ! Compute derivatives of (rho*u)
        call derivatives_ab(rho,rhox,rhoy, u,ux,uy,   a,ax,ay) !a=(rho*u)
        
        ! Compute derivatives of (rho*v)
        call derivatives_ab(rho,rhox,rhoy, v,vx,vy,   b,bx,by) !b=(rho*v)
        
        !-----------------------------------------------------------------------------
        
        ! Compute derivatives of (a*u)=(rho*u*u) !a=(rho*u)
        
        call derivatives_ab(a,ax,ay, u,ux,uy,      au,aux,auy)
        
        ! Compute derivatives of (a*v)=(rho*u*v) !a=(rho*u)
        
        call derivatives_ab(a,ax,ay, v,vx,vy,      av,avx,avy)
        
        ! Compute derivatives of (b*u)=(rho*v*u) !b=(rho*v)
        
        call derivatives_ab(b,bx,by,  u,ux,uy,     bu,bux,buy)
        
        ! Compute derivatives of (b*v)=(rho*v*v) !b=(rho*v)
        
        call derivatives_ab(b,bx, by, v,vx,vy,     bv,bvx,bvy)
        
        !-----------------------------------------------------------------------------
        
        ! Compute derivatives of (u*rH)
        
        call derivatives_ab( u,ux,uy, rhoH,rhoHx,rhoHy,  rhouH,rhouHx,rhouHy)
        
        ! Compute derivatives of (v*rH)
        
        call derivatives_ab( v,vx,vy, rhoH,rhoHx,rhoHy,  rhovH,rhovHx,rhovHy)
        
        !---------------------------------------------------------------------
        !---------------------------------------------------------------------

        !---------------------------------------------------------------------
        ! Store the inviscid terms in the forcing term array, f(:).
        !---------------------------------------------------------------------

        !------------------------------------------------------
        ! Continuity:         (rho*u)_x   +   (rho*v)_y
        S(1)  = (rhox*u + rho*ux) + (rhoy*v + rho*vy)

        !------------------------------------------------------
        ! Momentum:     (rho*u*u)_x + (rho*u*v)_y + px
        S(2)   =     aux     +    buy      + px
    
        !------------------------------------------------------
        ! Momentum:     (rho*u*v)_x + (rho*v*v)_y + px
        S(3)   =     avx     +    bvy      + py
    
        !------------------------------------------------------
        ! Momentum:     w is zero
        S(4)   = 0.0_p2
        !------------------------------------------------------
        ! Energy:       (rho*u*H)_x + (rho*v*H)
        S(5)  =    rhouHx   +   rhovHy
    
        if (iturb_type > TURB_INVISCID) then

            C0  = sutherland_constant/reference_temp
            xmr = M_inf/Re_inf
            ! mu  = xmr * (one + C0) / (Q(5) + C0)*Q(5)**(1.5_p2)
            
            ! txx = (2.0_p2/3.0_p2) * mu * (2.0_p2 * ux - vy ) ! wz = 0
            ! txy = mu * (uy + vx)

            ! tyx = txy
            ! tyy = (2.0_p2/3.0_p2) * mu * (2.0_p2 * vy - ux ) ! wz = 0


            ! ! Temperature
            ! call quotient_rule(p*gamma,px*gamma,py*gamma,rho,rhox,rhoy,T,dtdx,dtdy)
            ! dtdz = 0.0_p2

            ! ! 2nd Derivatives
            ! rxx = manufactured_sol(cr0,crs,crx,cry, 2,0,x,y)
            ! ryy = manufactured_sol(cr0,crs,crx,cry, 0,2,x,y)
            ! rxy = manufactured_sol(cr0,crs,crx,cry, 1,1,x,y)
            ! ryx = rxy

            ! uxx = manufactured_sol(cu0,cus,cux,cuy, 2,0,x,y)
            ! uyy = manufactured_sol(cu0,cus,cux,cuy, 0,2,x,y)
            ! uxy = manufactured_sol(cu0,cus,cux,cuy, 1,1,x,y)
            ! uyx = uxy

            ! vxx = manufactured_sol(cv0,cvs,cvx,cvy, 2,0,x,y)
            ! vyy = manufactured_sol(cv0,cvs,cvx,cvy, 0,2,x,y)
            ! vxy = manufactured_sol(cv0,cvs,cvx,cvy, 1,1,x,y)
            ! vyx = vxy
            
            ! pxx = manufactured_sol(cp0,cps,cpx,cpy, 2,0,x,y)
            ! pyy = manufactured_sol(cp0,cps,cpx,cpy, 0,2,x,y)
            ! pxy = manufactured_sol(cp0,cps,cpx,cpy, 1,1,x,y)
            ! pyx = pxy

            ! mudt = 0.5_p2 * mu * (1.0_p2 + 3.0_p2 * C0 / T) / (T + C0)
            ! mux = mudt * dtdx
            ! muy = mudt * dtdy
            ! muz = 0.0_p2

            ! ! d(pr)
            ! pxr = px * rho
            ! pyr = py * rho
            ! prx = p * rhox
            ! pry = p * rhoy
                
            ! ! dd(pr)
            ! pxxr = pxx * rho
            ! pxrx = px * rhox
            ! prxx = p  *  rxx

            ! pxyr = pxy * rho
            ! pxry = px * rhoy
            ! pyrx = py * rhox
            ! prxy = p *   rxy
            
            ! pyyr = pyy * rho
            ! pyry = py * rhoy
            ! pryy = p *   ryy

            ! pyxr = pxyr
            ! ! pyrx = py * rhox
            ! ! pxry = px * rhoy
            ! pryx = prxy

            ! ! Txx and Txy
            ! f  = gamma * (pxr - prx)
            ! fx = gamma * (pxxr - prxx)
            ! fy = gamma * (pxyr + pxry - pyrx - prxy)

            ! g = rho * rho
            ! gx = 2 * rho * rhox
            ! gy = 2 * rho * rhoy

            ! ! dummy is f/g which is Tx
            ! call quotient_rule(f,fx,fy,g,gx,gy,dummy,dtdxx,dtdxy)

            ! ! Tyx and Tyy
            ! f  = gamma * (pyr - pry)
            ! fx = gamma * (pyyr - pryy)
            ! fy = gamma * (pyxr + pyrx - pxry - pryx)

            ! ! g is the same
            ! call quotient_rule(f,fx,fy,g,gx,gy,dummy,dtdyx,dtdyy)

            ! qxx = ((-gamma)/(Pr*(gamma-1.0_p2))) * (mux * (-px/rho - p*rhox/rho**2) + &
            !      mu * (2.0_p2 * rhox*px / rho**2 + pxx/rho + p*rhox*rhox / rho**3 - p*rxx / rho**2) ) ! Eq 4.13.14 idlCFD
            ! qyy = ((-gamma)/(Pr*(gamma-1.0_p2))) * (muy * (-px/rho - p*rhox/rho**2) + &
            !      mu * (2.0_p2 * rhoy*py / rho**2 + pyy/rho + p*rhoy*rhoy / rho**3 - p*ryy / rho**2) ) ! Eq 4.13.15 idlCFD

            ! !d(tau)
            ! txx_x = (2.0_p2/3.0_p2) * ( mux * (2.0_p2 * ux - vy) + mu * (2.0_p2 * uxx - vyx) )
            ! txy_x = mux * (uy + vx) * mu * (uyx + vxx)
            ! tyx_x = txy_x

            ! tyy_y = (2.0_p2/3.0_p2) * ( muy * (2.0_p2 * vy - ux) + mu * (2.0_p2 * vyy - uxy) )
            ! txy_y = muy * (uy + vx) * mu * (uyy + vxy)
            ! tyx_y = txy_y

            ! tvx_x = u * txx_x + v * txy_x + txx * ux + txy * vx
            ! tvy_y = u * tyx_y + v * tyy_y + tyx * uy + tyy * vy

            ! F2hand = - txx_x - txy_y
            ! F3hand = - tyx_x - tyy_y
            ! F5hand = - tvx_x - tvy_y + qxx + qyy

            f2 = myF2(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            f3 = myF3(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            f5 = myF5(C0, Pr, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cu0, &
            cus, cux, cuy, cv0, cvs, cvx, cvy, gamma, gammamo, x, xmr, y)

            ! S(1) = S(1) ! no change to S(1)

            !------------------------------------------------------
            ! Momentum:     - (txx)_x - (txy)_y
            S(2) = S(2) + F2

            !------------------------------------------------------
            ! Momentum:     - (tyx)_x - (tyy)_y
            S(3) = S(3) + F3

            !------------------------------------------------------
            ! Momentum:    
            ! S(4) = S(4) ! 2D function => no change

            !------------------------------------------------------
            ! Energy:       (rho*u*H)_x + (rho*v*H)
            ! S(5)  = S(5) - tvx_x - tvy_y + qxx + qyy
            S(5)  = S(5) + F5

        endif
        !Return f.
    
        ! Note: Later, we'll perform the following to compute the residual for
        !                   dF(w)/dx + dG(w)/dy = f
        !       Step 1. Comptue the residual: Res=dF(w)/dx + dG(w)/dy.
        !       Step 2. Subtract f: Res = Res - f.
        !


    end subroutine fMMS

    !********************************************************************************
!* This function computes the sine function:
!*
!*       f =  a0 + as*sin(ax*x+ay*y)
!*
!* and its derivatives:
!*
!*     df/dx^nx/dy^ny = d^{nx+ny}(a0+as*sin(ax*x+ay*y))/(dx^nx*dy^ny)
!*
!* depending on the input parameters:
!*
!*
!* Input:
!*
!*  a0,as,ax,ay = coefficients in the function: f =  a0 + as*sin(ax*x+ay*y).
!*            x = x-coordinate at which the function/derivative is evaluated.
!*            y = y-coordinate at which the function/derivative is evaluated.
!*           nx = nx-th derivative with respect to x (nx >= 0).
!*           ny = ny-th derivative with respect to y (ny >= 0).
!*
!* Output: The function value.
!*
!*
!* Below are some examples:
!*
!*     f =  a0 + as*sin(ax*x+ay*y)            !<- (nx,ny)=(0,0)
!*
!*    fx =  ax * as*cos(ax*x+ay*y)            !<- (nx,ny)=(1,0)
!*    fy =  ay * as*cos(ax*x+ay*y)            !<- (nx,ny)=(0,1)
!*
!*   fxx = -ax**2 * as*sin(ax*x+ay*y)         !<- (nx,ny)=(2,0)
!*   fxy = -ax*ay * as*sin(ax*x+ay*y)         !<- (nx,ny)=(1,1)
!*   fyy = -ay**2 * as*sin(ax*x+ay*y)         !<- (nx,ny)=(0,2)
!*
!*  fxxx = -ax**3        * as*cos(ax*x+ay*y)  !<- (nx,ny)=(3,0)
!*  fxxy = -ax**2 *ay    * as*cos(ax*x+ay*y)  !<- (nx,ny)=(2,1)
!*  fxyy = -ax    *ay**2 * as*cos(ax*x+ay*y)  !<- (nx,ny)=(1,2)
!*  fyyy = -       ay**3 * as*cos(ax*x+ay*y)  !<- (nx,ny)=(0,3)
!*
!* fxxxx =  ax**4        * as*sin(ax*x+ay*y)  !<- (nx,ny)=(4,0)
!* fxxxy =  ax**3 *ay    * as*sin(ax*x+ay*y)  !<- (nx,ny)=(3,1)
!* fxxyy =  ax**2 *ay**2 * as*sin(ax*x+ay*y)  !<- (nx,ny)=(2,2)
!* fxyyy =  ax    *ay**3 * as*sin(ax*x+ay*y)  !<- (nx,ny)=(1,3)
!* fyyyy =         ay**4 * as*sin(ax*x+ay*y)  !<- (nx,ny)=(0,4)
!*
!* and so on.
!*
!*
!********************************************************************************
 function manufactured_sol(a0,as,ax,ay, nx,ny,x,y) result(fval)

    implicit none
   
   !Input
    real(p2), intent(in) :: a0, as, ax, ay, x, y
    integer , intent(in) :: nx, ny
   
   !Output
    real(p2)             :: fval
   
     if (nx < 0 .or. ny < 0) then
      write(*,*) " Invalid input: nx and ny must be greater or equal to zero... Try again."
      stop
     endif
   
     if ( nx+ny == 0 ) then
   
      fval = a0 + as*sin(ax*x + ay*y)
   
     elseif ( mod(nx+ny,2) == 0 ) then
   
      fval = - (ax**nx * ay**ny)*as*sin(ax*x + ay*y)
      if ( mod(nx+ny,4)   == 0 ) fval = -fval
   
     else
   
      fval = (ax**nx * ay**ny)*as*cos(ax*x + ay*y)
      if ( mod(nx+ny+1,4) == 0 ) fval = -fval
   
     endif
   
   
    end function manufactured_sol
    
    !********************************************************************************

    !********************************************************************************
    !
    ! This subroutine computes first derivatives of a quadratic term
    !
    !  Input: a, ax, ay !Function value a, and its derivatives, (ax,ay).
    !         b, bx, by !Function value b, and its derivatives, (bx,by).
    !
    ! Output: ab = a*b, abx = d(a*b)/dx, aby = d(a*b)/dy.
    !
    !********************************************************************************
    subroutine derivatives_ab(a,ax,ay,  b,bx,by, ab,abx,aby)

        implicit none
    
        !Input
        real(p2), intent( in) ::  a,  ax,  ay
        real(p2), intent( in) ::  b,  bx,  by
    
        !Output
        real(p2), intent(out) :: ab, abx, aby
    
        ab    = a*b 
        abx   = ax*b + a*bx
        aby   = ay*b + a*by
    
    end subroutine derivatives_ab

    ! Compute the quotient rule h' = (f'g - fg')/g^2
    subroutine quotient_rule(f, fx, fy, g, gx, gy, h, hx, hy)
        real(p2), intent(in) :: f, fx, fy
        real(p2), intent(in) :: g, gx, gy

        real(p2), intent(out):: h, hx, hy

        h = f/g
        hx = (fx*g - f*gx) / (g**2)
        hy = (fy*g - f*gy) / (g**2)

    end subroutine

end module mms