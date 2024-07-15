module limiter

    implicit none

    public :: compute_limiter

    contains

    subroutine compute_limiter

        use common          , only : p2, zero
        
        use grid            , only : ncells, cell, x, y, z, cell
        
        use solution        , only : ccgradq, phi, q
      
        ! use least_squares   , only : lsq 

        implicit none
        ! Some local vars
        integer  :: i, ivar, k, nghbr_cell, iv
        real(p2) :: qmin, qmax, xc, yc, zc, xp, yp, zp, qf, dqm, dqp
        real(p2) :: phi_vertex, phi_vertex_min, limiter_beps
        real(p2) :: phi_var_min

        !allocate(phi(ncells)) ! possible memory leak? Moved allocation to steady solve subroutine (only called once)
        limiter_beps = 1.0e-14_p2
        !loop over cells
        cell_loop : do i = 1,ncells
            variable_loop : do ivar = 1,5
                qmin = q(ivar,i)
                qmax = q(ivar,i)
                nghbr_loop : do k = 1,cell(i)%nnghbrs
                    nghbr_cell = cell(i)%nghbr(k)
                    qmin = min(qmin, q(ivar,nghbr_cell) )
                    qmax = max(qmax, q(ivar,nghbr_cell) )
                end do nghbr_loop
                ! Compute phi to enforce maximum principle at vertices (MLP)
                xc = cell(i)%xc
                yc = cell(i)%yc
                zc = cell(i)%zc

                ! Loop over vertices of the cell
                vertex_loop : do k = 1,cell(i)%nvtx
                    iv = cell(i)%vtx(k)
                    xp = x(iv)
                    yp = y(iv)
                    zp = z(iv)

                    ! Linear reconstruction to the vertex k
                    qf = q(ivar,i) + ccgradq(1,ivar,i)*(xp-xc) + &
                                     ccgradq(2,ivar,i)*(yp-yc) + &
                                     ccgradq(3,ivar,i)*(zp-zc)

                    ! compute dq^-
                    dqm = qf - q(ivar,i)

                    !Compute dq^+.
                    if ( dqm > zero ) then
                        dqp = qmax - q(ivar,i)
                    else
                        dqp = qmin - q(ivar,i)
                    endif

                    ! Limiter function: Venkat limiter

                    phi_vertex = vk_limiter(dqp, dqm, cell(i)%vol)
 
                    ! Keep the minimum over the control points (vertices).
                    if (k==1) then
                        phi_vertex_min = phi_vertex
                    else
                        phi_vertex_min = min(phi_vertex_min, phi_vertex)
                    endif
                end do vertex_loop
                if (ivar == 1) then
                    phi_var_min = phi_vertex_min
                else
                    phi_var_min = min(phi_var_min, phi_vertex_min)
                endif
            end do variable_loop
            phi(i) = phi_var_min
        end do cell_loop
                    
    end subroutine compute_limiter

    !********************************************************************************
    !* -- Venkat Limiter Function--
    !*
    !* 'Convergence to Steady State Solutions of the Euler Equations on Unstructured
    !*  Grids with Limiters', V. Venkatakrishnan, JCP 118, 120-130, 1995.
    !*
    !* The limiter has been implemented in such a way that the difference, b, is
    !* limited in the form: b -> vk_limiter * b.
    !*
    !* ------------------------------------------------------------------------------
    !*  Input:     a, b     : two differences
    !*
    !* Output:   vk_limiter : to be used as b -> vk_limiter * b.
    !* ------------------------------------------------------------------------------
    !*
    !
    !
    ! Note: This is unaltered from the edu_euler code.  We'll see if I have any need to edit it in the future
    !********************************************************************************
    pure function vk_limiter(a, b, vol)

        use common    , only : p2, two, pi, half, six, third

        real(p2), intent(in) :: a, b
        real(p2), intent(in) :: vol

        real(p2)             :: vk_limiter
        real(p2)             :: Kp
        real(p2)             :: eps2, diameter

        !  real(p2)             :: r

        Kp = 5.0_p2   !<<<<< Adjustable parameter K

        !   Mesh dependent constant (See Eqn.(33) in the paper) in Venkat limiter.
        !   chokkei = (6.0_p2*elm(i)%vol/pi)**(1.0_p2/3.0_p2) ! 3D version
        !      eps2 = (Kp*chokkei)**3

        diameter = six*(vol/pi)**third  ! 2D version = 2 times the diamater
            eps2 = (Kp*diameter)**3

        ! This is the form used by Venkat. This is in the form of
        ! limited_slope(a,b) / b, so that limited_slope(a,b) part resembles
        ! Van Albada's original limiter. And this way, he follows Van Albada
        ! and introduced epsilon to avoid limiting in nearly constant regions.
        !
        !    vk_limiter = ( b*(a**2 + eps2) + two*b**2*a )/(a**2 + two*b**2 + a*b + eps2) / b

        ! The above is equivalent to the following. This is within [0,1], well,
        ! it overshoots 1.0 near r=1.0, but approaches to 1.0 as r goes large...

        vk_limiter = ( (a**2 + eps2) + two*b*a )/(a**2 + two*b**2 + a*b + eps2)

        !            r = a/b
        !   vk_limiter = ( r + abs(r) ) / (one + abs(r) )
        !   vk_limiter = ( two*r ) / (one + r )
        !   vk_limiter = two*( a*b + eps )/( a**2 + b**2 + two*eps ) / b

  end function vk_limiter

end module limiter