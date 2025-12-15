module mms

    use common , only : p2, pi

    implicit none

    private
    public fMMS, txx_x, tyy_y, txy_x, txy_y

    ! coefficients
    real(p2), parameter :: arx = 1.5_p2
    real(p2), parameter :: d = 0.25
    real(p2), parameter :: pi2 = 2.0_p2 * pi

    real(p2) :: txx_x, tyy_y, txy_x, txy_y

    contains

    subroutine fMMS(x,y,z, Q,S, gradQ)

        use config , only : M_inf, Pr, sutherland_constant, Re_inf, M_inf, reference_temp, mms_include

        use utils , only : iturb_type, TURB_INVISCID

        use common , only : half, one

        use solution , only : w2q, gamma, gammamo

        use mms_funcs

        implicit none

        real(p2),                 intent(in)  :: x,y,z ! coordinates
        real(p2), dimension(5),   intent(out) :: Q, S  ! Solution and source terms
        real(p2), dimension(3,5), optional, intent(inout) :: gradQ

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

        real(p2) :: C0, xmr, a2

        real(p2), dimension(5) :: wtmp

        S(:) = 0.0_p2 ! avoid undefined behavior if mms_include(1) = .false.
        
       !-----------------------------------------------------------
        ! Constants for the exact solution: c0 + cs*sin(cx*x+cy*y).
        !
        ! Note: Make sure the density and pressure are positive (c0 > cs).
        ! Note: These values are passed to the subroutine:
        !         manufactured_sol(c0,cs,cx,cy, nx,ny,x,y),
        !       whcih returns the solution value or derivatives.

        !-----------------------------------------
        ! Density    = cr0 + crs*sin(crx*x+cry*y)

        cr0 =  1.12_p2
        crs =  0.15_p2
        crx =  3.12_p2*pi !/ 1024.0_p2
        cry =  2.92_p2*pi !/ 1024.0_p2
        ! crx =  0.000001_p2
        ! cry =  0.000001_p2
        ! crs = 1.0_p2 / crx

        !-----------------------------------------
        ! X-velocity = cu0 + cus*sin(cux*x+cuy*y)

        cu0 =  1.32_p2
        cus =  0.06_p2
        cux =  2.09_p2*pi !/ 1024.0_p2
        cuy =  3.12_p2*pi !/ 1024.0_p2
        ! cux =  0.000001_p2
        ! cuy =  0.000001_p2
        ! cus = 1.0_p2 / cux

        !-----------------------------------------
        ! Y-velocity = cv0 + cvs*sin(cvx*x+cvy*y)

        cv0 =  1.18_p2
        cvs =  0.03_p2
        cvx =  2.15_p2*pi !/ 1024.0_p2
        cvy =  3.32_p2*pi !/ 1024.0_p2
        ! cvx =  0.000001_p2
        ! cvy =  0.000001_p2
        ! cvs = 1.0_p2 / cvx

        !-----------------------------------------
        ! Pressure   = cp0 + cps*sin(cpx*x+cpy*y)

        cp0 =  1.62_p2
        cps =  0.31_p2
        cpx =  3.79_p2*pi !/ 1024.0_p2
        cpy =  2.98_p2*pi !/ 1024.0_p2
        ! cpx =  0.000001_p2
        ! cpy =  0.000001_p2
        ! cps = 1.0_p2 / cpx

        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------
        ! Part I: Compute w = [rho,u,v,p] and grad(w).
        !-----------------------------------------------------------------------------
        !-----------------------------------------------------------------------------

        !------------------------------------------------------------------------
        ! rho: Density and its 1st derivatives

        rho  = manufactured_sol(cr0,crs,crx,cry, 0,0,x,y)
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

        if (present(gradQ)) then
            gradQ = 0.0_p2
            gradQ(:,1) = (/px, py, 0.0_p2/)
            gradQ(:,2) = (/ux, uy, 0.0_p2/)
            gradQ(:,3) = (/vx, vy, 0.0_p2/)
            gradQ(:,4) = 0.0_p2
            
            a2 = gamma*p/rho
            
            ! gradQ(:,5)    = ( gamma*gradQ(:,1) - a2*(/rhox, rhoy, 0.0_p2/) )/rho
            gradQ(:,5) = (/dTx(cp0, cps, cpx, cpy, cr0, crs, crx, cry, gamma, x, y), &
                          dTy(cp0, cps, cpx, cpy, cr0, crs, crx, cry, gamma, x, y), &
                          0.0_p2/)
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

        if (mms_include(1)) then
            ! ------------------------------------------------------
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
        endif
            
        if ((iturb_type > TURB_INVISCID) .and. (mms_include(2).or.mms_include(3)) ) then

            C0  = sutherland_constant/reference_temp
            xmr = M_inf/Re_inf

            txx_x = mytxx_x(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            tyy_y = my_tyy_y(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            txy_x  = my_txy_x(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            txy_y = my_txy_y(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            F2 = myF2(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            ! F3 = myF3(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            ! cuy, cvs, cvx, cvy, gamma, x, xmr, y)
            
            F3 = myNewF3(C0, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cus, cux, &
            cuy, cvs, cvx, cvy, gamma, x, xmr, y)

            ! F5 = myF5(C0, Pr, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cu0, &
            ! cus, cux, cuy, cv0, cvs, cvx, cvy, gamma, gammamo, x, xmr, y)

            F5 = mynewF5(C0, Pr, cp0, cps, cpx, cpy, cr0, crs, crx, cry, cu0, &
            cus, cux, cuy, cv0, cvs, cvx, cvy, gamma, gammamo, x, xmr, y)

            ! S(1) = S(1) ! no change to S(1)

            !------------------------------------------------------
            ! Momentum:     - (txx)_x - (txy)_y
            S(2) = S(2) + F2
            ! S(2) = S(2) - txx_x
            ! S(2) = S(2) - txy_y

            !------------------------------------------------------
            ! Momentum:     - (tyx)_x - (tyy)_y
            S(3) = S(3) + F3
            ! S(3) = S(3) - tyy_y
            ! S(3) = S(3) - txy_x

            !------------------------------------------------------
            ! Momentum:    
            ! S(4) = S(4) ! 2D function => no change

            !------------------------------------------------------
            ! Energy:       (rho*u*H)_x + (rho*v*H)
            ! S(5)  = S(5) - tvx_x - tvy_y + qxx + qyy
            S(5)  = S(5) + F5
            ! S(5)  = S(5) + my_dq(C0, Pr, cp0, cps, cpx, cpy, cr0, crs, crx, cry, gamma, &
                                ! gammamo, x, xmr, y)

        endif
        !Return f.
    
        ! Note: Later, we'll perform the following to compute the residual for
        !                   dF(w)/dx + dG(w)/dy = f
        !       Step 1. Comptue the residual: Res=dF(w)/dx + dG(w)/dy.
        !       Step 2. Subtract f: Res = Res - f.
        !

        
    end subroutine fMMS

    function manufactured_sol(a0,as,ax,ay, nx,ny,x,y) result(fval)
        
        use utils , only : imms_type, MMS_LIN, MMS_QUAD, MMS_SIN 
        implicit none

        !Input
        real(p2), intent(in) :: a0, as, ax, ay, x, y
        integer , intent(in) :: nx, ny

        !Output
        real(p2)             :: fval

        select case(imms_type)
        case(MMS_LIN)
            fval = manufactured_sol_linear(a0,as,ax,ay, nx,ny,x,y)
        case(MMS_QUAD)
            fval = manufactured_sol_quadratic(a0,as,ax,ay, nx,ny,x,y)
        case(MMS_SIN)
            fval = manufactured_sol_sin(a0,as,ax,ay, nx,ny,x,y)
        case default
            write(*,*) 'imms_type: ', imms_type , 'invalid. Stopping.'
            stop
        end select

    end function manufactured_sol

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
    function manufactured_sol_sin(a0,as,ax,ay, nx,ny,x,y) result(fval)

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
   
   
    end function manufactured_sol_sin

     !********************************************************************************
    !* This function computes the linear function (this should be exact):
    !*
    !*       f =  a0 + ax*x+ay*y (as is not used)
    !*
    !* and its derivatives:
    !*
    !*     df/dx  = ax
    !*     df/dy  = ay
    !*     d2/dn2 = 0
    !*
    !* depending on the input parameters:
    !*
    !*
    !* Input:
    !*
    !*     a0,ax,ay = coefficients in the function: f =  a0 + as*sin(ax*x+ay*y).
    !*            x = x-coordinate at which the function/derivative is evaluated.
    !*            y = y-coordinate at which the function/derivative is evaluated.
    !*           nx = nx-th derivative with respect to x (nx >= 0).
    !*           ny = ny-th derivative with respect to y (ny >= 0).
    !*
    !* Output: The function value.
    !*
    !* Note: While the functions for P, Rho, U, and V are linear, the resulting function for 
    !*       T will not be.  Instead T = P*gamma/Rho.  As a result, the gradient for T will 
    !*       not be exact.
    !*
    !*
    !********************************************************************************
   
    function manufactured_sol_linear(a0,as,ax,ay, nx,ny,x,y) result(fval)

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
   
      fval = a0 + ax*x + ay*y
   
     elseif ( nx+ny > 1 ) then
   
      fval = 0
      
     elseif (nx == 1) then
   
      fval = ax

     else
      fval = ay
     endif
   
   
    end function manufactured_sol_linear
    

    !********************************************************************************
    !* This function computes the quadratic function:
    !*
    !*       f =  a0 + ax*x*x+ay*y*y + as*x*y (as is not used)
    !*
    !* and its derivatives:
    !*
    !*     df/dx  = 2*ax*x + as*y
    !*     df/dy  = 2*ay*y + as*x
    !*     df/dx2 = 2ax
    !*     df/dy2 = 2ay
    !*     df/dxy = as
    !*
    !* depending on the input parameters:
    !*
    !*
    !* Input:
    !*
    !*     a0,ax,ay = coefficients in the function: f =  a0 + as*sin(ax*x+ay*y).
    !*            x = x-coordinate at which the function/derivative is evaluated.
    !*            y = y-coordinate at which the function/derivative is evaluated.
    !*           nx = nx-th derivative with respect to x (nx >= 0).
    !*           ny = ny-th derivative with respect to y (ny >= 0).
    !*
    !* Output: The function value.
    !*
    !*
    !********************************************************************************
   
    function manufactured_sol_quadratic(a0,as,ax,ay, nx,ny,x,y) result(fval)

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
   
      fval = a0 + ax*x*x + ay*y*y + as*x*y
   
     elseif ( nx+ny > 2 ) then
   
      fval = 0
     
     elseif (nx == 1 .and. ny == 0) then
   
      fval = 2.0_p2 * ax * x + as * y

     elseif (nx == 0 .and. ny == 1) then
   
      fval = 2.0_p2 * ay * y + as * x

     elseif (nx == 1 .and. ny == 1) then
   
      fval = as

     elseif (nx == 2) then
   
      fval = 2.0_p2 * ax

     else ! (nx == 2) then
   
      fval = 2.0_p2 * ay

     endif
   
   
    end function manufactured_sol_quadratic
    
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
        
        implicit none
        
        real(p2), intent(in) :: f, fx, fy
        real(p2), intent(in) :: g, gx, gy

        real(p2), intent(out):: h, hx, hy

        h = f/g
        hx = (fx*g - f*gx) / (g**2)
        hy = (fy*g - f*gy) / (g**2)

    end subroutine

end module mms