module res_sa

    use common , only : p2

    use sa_vars
    
    implicit none

    private

    public compute_res_sa



    contains

    subroutine compute_res_sa

        use common  , only : zero, half, one

        use config  , only : use_limiter, CFL_turb

        use utils   , only : ibc_type

        use grid    , only : ncells, cell,  &
                             nfaces, face,  &
                             nb,     bound, &
                             face_nrml,     &
                             face_nrml_mag, &
                             face_centroid
        
        use gradient , only : compute_gradient_turb

        use turb     , only : turb_res, turb_jac, turb_var, phi_turb, ccgrad_turb_var, vgrad_turb_var

        use solution_vars , only : ccgradq, q, kth_nghbr_of_1, kth_nghbr_of_2, wsn

        use solution , only : q2u, q2rho

        use turb_bc , only : sa_rhstate

        use bc_states, only : get_right_state

        use limiter  , only : compute_limiter_turb

        use wall_distance , only : cell_wall_distance

        use direct_solve , only : safe_invert_scalar

        use viscosity , only : compute_viscosity

        ! Grid Vars
        integer                     :: cell1, cell2
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid
        real(p2)                    :: face_mag
        real(p2)                    :: nut1, nut2, nutb
        real(p2), dimension(5)      :: q1, q2, qb
        real(p2), dimension(3,5)    :: gradq
        real(p2), dimension(3)      :: gradnut1, gradnut2, gradnutb
        real(p2)                    :: phi1, phi2
        real(p2)                    :: d1
        real(p2)                    :: nu1
        real(p2)                    :: rho1
        real(p2)                    :: itwsn, dtaui

        ! real(p2)                    :: num_flux, num_jac1, num_jac2
        real(p2)                    :: num_jacsrc
        real(p2), dimension(2)      :: num_flux, num_jac1, num_jac2 ! convective terms are cell specific
        real(p2)                    :: nsource

        integer :: iface, ib, icell
        integer :: k, nk
        integer :: face_sides

        turb_res(:,:) = zero
        
        do icell = 1,ncells
            turb_jac(icell,1)%diag = zero
            turb_jac(icell,1)%off_diag(:) = zero
        end do

        ! update turbulent variable gradients
        call compute_gradient_turb(0)

        ! hmmm I guess I still gotta sort the limiter out
        call compute_limiter_turb

        loop_faces : do iface = 1,nfaces
            cell1 = face(1,iface)
            cell2 = face(2,iface)
            q1 = q(1:5, cell1)
            q2 = q(1:5, cell2)
            nut1 = turb_var(cell1,1)
            nut2 = turb_var(cell2,1)
            gradnut1 = ccgrad_turb_var(:,cell1,1)
            gradnut2 = ccgrad_turb_var(:,cell2,1)

            ! Face normal
            unit_face_normal = face_nrml(1:3,iface)

            ! Limiters
            if (use_limiter) then
                phi1 = phi_turb(cell1)
                phi2 = phi_turb(cell2)
            else 
                phi1 = one
                phi2 = one
            end if

            call sa_invFlux(                    nut1,     nut2, &
                                                  q1,       q2, &
                                            gradnut1, gradnut2, &
                                              unit_face_normal, &
                cell(cell1)%xc, cell(cell1)%yc, cell(cell1)%zc, & !<- Left  cell centroid
                cell(cell2)%xc, cell(cell2)%yc, cell(cell2)%zc, & !<- Right cell centroid
                                        face_centroid(1,iface), &
                                        face_centroid(2,iface), &
                                        face_centroid(3,iface), & !<- face midpoint
                                             phi1,        phi2, & !<- Limiter functions
                               num_flux, num_jac1, num_jac2  ) !<- Output

            !Cell 1
            turb_res(cell1,1) =             turb_res(cell1,1)             + num_flux(1) * face_nrml_mag(iface)

            turb_jac(cell1,1)%diag =        turb_jac(cell1,1)%diag        + num_jac1(1) * face_nrml_mag(iface)
            k = kth_nghbr_of_1(iface)
            turb_jac(cell1,1)%off_diag(k) = turb_jac(cell1,1)%off_diag(k) + num_jac1(2) * face_nrml_mag(iface)
            
            ! Cell 2 Note the + sign because the flux subroutine already accounts sign
            turb_res(cell2,1) =             turb_res(cell2,1)             + num_flux(2) * face_nrml_mag(iface)

            turb_jac(cell2,1)%diag =        turb_jac(cell2,1)%diag        + num_jac2(1) * face_nrml_mag(iface)
            k = kth_nghbr_of_2(iface)
            turb_jac(cell2,1)%off_diag(k) = turb_jac(cell2,1)%off_diag(k) + num_jac2(2) * face_nrml_mag(iface)

            ! Add contribution to the wave speed which (assuming I understand the Rankine Hugonot relation correctly)
            ! is equivalent to the jacobian of the convective term.
            ! One will always be zero so we take the absolute value of whichever one isn't.
            ! itwsn = max(abs(num_jac1),abs(num_jac2))

            ! Diffusion Flux terms
            call sa_viscFlux(                   nut1,     nut2, &
                                                  q1,       q2, &
                                            gradnut1, gradnut2, &
                                              unit_face_normal, &
                cell(cell1)%xc, cell(cell1)%yc, cell(cell1)%zc, & !<- Left  cell centroid
                cell(cell2)%xc, cell(cell2)%yc, cell(cell2)%zc, & !<- Right cell centroid
                                  num_flux, num_jac1, num_jac2  ) !<- Output

            !Cell 1
            turb_res(cell1,1)             = turb_res(cell1,1)             + num_flux(1) * face_nrml_mag(iface)

            turb_jac(cell1,1)%diag        = turb_jac(cell1,1)%diag        + num_jac1(1) * face_nrml_mag(iface)
            k = kth_nghbr_of_1(iface)
            turb_jac(cell1,1)%off_diag(k) = turb_jac(cell1,1)%off_diag(k) + num_jac1(2) * face_nrml_mag(iface)
            
            ! Cell 2
            turb_res(cell2,1)             = turb_res(cell2,1)             + num_flux(2) * face_nrml_mag(iface)

            turb_jac(cell2,1)%diag        = turb_jac(cell2,1)%diag        + num_jac2(1) * face_nrml_mag(iface)
            k = kth_nghbr_of_2(iface)
            turb_jac(cell2,1)%off_diag(k) = turb_jac(cell2,1)%off_diag(k) + num_jac1(2) * face_nrml_mag(iface)
                      
        end do loop_faces
        
        gradnut2 = zero
        
        bound_loop : do ib = 1,nb
            bfaces_loop : do iface = 1,bound(ib)%nbfaces
                cell1 = bound(ib)%bcell(iface)
                
                gradnut1 = ccgrad_turb_var(:,cell1,1)
                
                bface_centroid   = bound(ib)%bface_center(:,iface)
                unit_face_normal = bound(ib)%bface_nrml(:,iface)
                face_mag         = bound(ib)%bface_nrml_mag(iface)

                nut1 = turb_var(cell1,1)
                q1   =        q(:,cell1)

                if (use_limiter) then
                    phi1 = phi_turb(cell1)
                    phi2 = one
                else 
                    phi1 = one
                    phi2 = one
                end if

                call sa_rhstate(nut1, ibc_type(ib), nutb)

                call get_right_state(q1,unit_face_normal,ibc_type(ib),qb)

                call sa_invFlux(                    nut1,     nutb, &
                                                      q1,       qb, &
                                                gradnut1, gradnut2, &
                                                  unit_face_normal, &
                    cell(cell1)%xc, cell(cell1)%yc, cell(cell1)%zc, & !<- Left  cell centroid
             bface_centroid(1),bface_centroid(2),bface_centroid(3), & !<- Face midpoint
             bface_centroid(1),bface_centroid(2),bface_centroid(3), & !<- Face midpoint
                                                 phi1,        phi2, & !<- Limiter functions
                                      num_flux, num_jac1, num_jac2  ) !<- Output

                !Cell 1 only
                turb_res(cell1,1)      = turb_res(cell1,1)             + num_flux(1) * face_mag

                turb_jac(cell1,1)%diag = turb_jac(cell1,1)%diag        + num_jac1(1) * face_mag
                ! No off diagonal terms and the second term of num flux is ignored.
                
                face_sides = bound(ib)%bfaces(1,iface)

                gradnutb = zero
                do k = 1,face_sides
                    nk = bound(ib)%bfaces(k+1,iface)
                    gradnutb = gradnutb + vgrad_turb_var(:,nk,1)
                end do
                gradnutb = gradnutb / real(face_sides,p2)

                ! Diffusion Flux terms
                call sa_viscFlux(                   nut1,     nutb, &
                                                      q1,       qb, &
                                                gradnutb, gradnutb, & !We want gradface = gradnutb
                                                  unit_face_normal, &
                    cell(cell1)%xc, cell(cell1)%yc, cell(cell1)%zc, & !<- Left  cell centroid
             bface_centroid(1),bface_centroid(2),bface_centroid(3), & !<- Face midpoint
                                      num_flux, num_jac1, num_jac2  ) !<- Output

                !Cell 1
                turb_res(cell1,1)      = turb_res(cell1,1)             + num_flux(1) * face_mag

                turb_jac(cell1,1)%diag = turb_jac(cell1,1)%diag        + num_jac1(1) * face_mag

            end do bfaces_loop

        end do bound_loop

        ! Source loop
        do icell     = 1,ncells
            nut1     = turb_var(icell,1)
            d1       = cell_wall_distance(icell)
            gradq    = ccgradq(:,:,icell)
            gradnut1 = ccgrad_turb_var(:,icell,1)
            rho1 = q2rho(q(:,icell))
            nu1 = compute_viscosity(q(5,icell)) / rho1 ! kinematic viscosity
            call sa_source(nut1,d1,nu1, gradq, gradnut1, nsource, num_jacsrc)

            ! We have to subtract the source term to move it to the LHS
            turb_res(icell,1)      = turb_res(icell,1) - nsource * cell(icell)%vol

            turb_jac(icell,1)%diag = turb_jac(icell,1)%diag - num_jacsrc * cell(icell)%vol

        end do

        do icell = 1,ncells

            ! TODO add pseudo-transient term to this
            dtaui = CFL_turb * cell(icell)%vol/( half * wsn(icell) )
            turb_jac(icell,1)%diag = turb_jac(icell,1)%diag + cell(icell)%vol / dtaui
            ! turb_jac(icell,1)%diag = turb_jac(icell,1)%diag + half * wsn(icell) / CFL_turb

            turb_jac(icell,1)%diag_inv = safe_invert_scalar(turb_jac(icell,1)%diag)
        end do

    end subroutine compute_res_sa

    subroutine sa_invFlux(nut1, nut2, q1, q2, gradnut1, gradnut2, n12, xc1, yc1, zc1, xc2, yc2, zc2, &
                             xm, ym, zm, phi1, phi2, nut_flux, jac1, jac2 )

        use config , only : rans_accuracy

        use common , only : half, zero

        use solution_vars , only : nq, iu, iv, iw
                            
        implicit none

        real(p2),               intent(in) :: nut1, nut2
        real(p2), dimension(:), intent(in) :: q1,  q2
        real(p2), dimension(:), intent(in) :: gradnut1, gradnut2
        real(p2), dimension(:), intent(in) :: n12
        real(p2),               intent(in) :: xc1, yc1, zc1, xc2, yc2, zc2
        real(p2),               intent(in) :: xm, ym, zm
        real(p2),               intent(in) :: phi1, phi2

        real(p2), dimension(2), intent(out):: nut_flux
        real(p2), dimension(2), intent(out):: jac1, jac2

        real(p2) :: nutL, nutR

        real(p2), dimension(nq) :: qL, qR, qi
        real(p2)                :: vP, vM, vBar

        real(p2), parameter :: eig_min = 1e-06


        if (rans_accuracy == 2) then
            nutL = nut1 + phi1 * ( gradnut1(1)*(xm-xc1) + gradnut1(2)*(ym-yc1) + gradnut1(3)*(zm-zc1) ) ! gradnut <=> gradnut (var) 
            nutR = nut2 + phi2 * ( gradnut2(1)*(xm-xc2) + gradnut2(2)*(ym-yc2) + gradnut2(3)*(zm-zc2) ) ! u <=> nut (vars)
        else
            nutL = nut1 
            nutR = nut2
        end if

        qL = q1
        qR = q2

        vBar = qL(iu) * n12(1) + qL(iv) * n12(2) + qL(iw) * n12(3)
        vP   = half * (vBar + abs(vBar) )
        vM   = half * (vBar - abs(vBar) )
        
        nut_flux(1) = vP * nutL + vM * nutR
        jac1(1)       = vP ! diag (dF/dnuL)
        jac1(2)       = vM ! off diag (dF/dnuR)

        vBar = -( qR(iu) * n12(1) + qR(iv) * n12(2) + qR(iw) * n12(3) )
        vP   = half * (vBar + abs(vBar) )
        vM   = half * (vBar - abs(vBar) )
        
        nut_flux(2) = vP * nutR + vM * nutL
        jac2(1)     = vP ! diag (dF/dnuR)
        jac2(2)     = vM ! off diag (dF/dnuL)
        ! Apply Harten's Entropy fix as described by Eq 6 in: https://doi.org/10.2514/1.J058549
        ! if (abs(vF) < eig_min) then
        !     vF = half * ( (vF*vF / eig_min) + eig_min) * sign(one,vF)
        ! endif

    end subroutine sa_invFlux

    subroutine sa_viscFlux(nut1,nut2,q1,q2,gradnut1,gradnut2,n12,xc1,yc1,zc1,xc2,yc2,zc2, nut_flux, jac1, jac2 )

        use common , only : half

        use solution_vars,only : ndim, nq

        use solution , only : q2u

        use viscosity,only: compute_viscosity
        implicit none

        real(p2),               intent(in) :: nut1, nut2
        real(p2), dimension(:), intent(in) :: gradnut1, gradnut2
        real(p2), dimension(:), intent(in) :: n12
        real(p2),               intent(in) :: xc1, yc1, zc1, xc2, yc2, zc2
        real(p2), dimension(:), intent(in) :: q1,q2
        
        real(p2), dimension(2), intent(out):: nut_flux
        real(p2), dimension(2), intent(out):: jac1, jac2

        ! Local
        real(p2), dimension(ndim)   :: gradnut_face
        real(p2), dimension(ndim)   :: ds,  dsds2
        real(p2)                    :: T, rho
        real(p2), dimension(nq)     :: u
        real(p2)                    :: muf ! dynamic viscosity
        real(p2)                    :: nuf ! kinematic viscosity
        real(p2)                    :: nutf ! face nut
        real(p2)                    :: normal_face_grad
        real(p2)                    :: term1, term21, term22

        ! Calculate the face gradients
        ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
        dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds(:)/ds**2

        gradnut_face(:) = half * (gradnut1(:) + gradnut2(:))
        gradnut_face(:) = gradnut_face(:) + ( (nut2 - nut1) - dot_product(gradnut_face(:),ds)) * dsds2
        
        T   = half * ( q1( 5 ) + q2( 5 ) )
        u   = half * ( q2u(q1) + q2u(q2) )
        nutf = half * ( nut1    + nut2    )
        rho = u(1)
        muf  = compute_viscosity(T)
        nuf  = muf / rho ! Kinematic Viscosity

        term1 = (one + cb2) * (nuf + nutf)
        term21 = cb2 * (nuf + nut1) ! cross diffusion term doesn't use the face gradient it uses the cell value
        term22 = cb2 * (nuf + nut2)

        normal_face_grad = dot_product( gradnut_face, n12 )

        nut_flux(1) = iSIGMA * (term1 - term21) * normal_face_grad
        jac1(:)     = iSIGMA * (one + cb2) * normal_face_grad * half
        jac1(1)     = (jac1(1) - cb2 * iSIGMA * normal_face_grad)

        nut_flux(2) = -( iSIGMA * (term1 - term22) * normal_face_grad )
        jac2(:)     = - ( iSIGMA * (one + cb2) * normal_face_grad ) * half
        jac2(1)     =   ( jac2(1) + cb2 * iSIGMA * normal_face_grad ) ! have to be a little careful with the signs here

    end subroutine sa_viscFlux

    subroutine sa_source(nut,distance,kvisc,gradQ,gradnut, source,dsource)
        
        use common , only : zero, sixth, fivesixth, one, three, six, & 
                            ix, iy, iz
        
        use solution_vars , only : iu, iv, iw

        implicit none

        real(p2),                 intent(in)  :: nut, distance, kvisc
        real(p2), dimension(:,:), intent(in)  :: gradQ
        real(p2), dimension(:),   intent(in)  :: gradnut

        real(p2),                 intent(out) :: source, dsource

        ! SA Variables
        real(p2) :: Chi
        real(p2) :: fv1, fv2, fw, ft2
        real(p2) :: g, r, Omega, Shat, Sbar
        ! real(p2) :: Wij
        real(p2) :: uy, uz, vx, vz, wx, wy

        ! SA Jacobian Variables
        real(p2) :: dChi
        real(p2) :: dfv1, dfv2, dfw, dft2
        real(p2) :: dg, dr, dShat, dSbar

        ! Source terms
        real(p2) :: prod, dest, s1
        real(p2) :: dprod, ddest
        ! Terms from Eq 15 and 16 of Spalart 1992 and their derivative
        real(p2) :: p, d, pprm, dprm

        ! Temp vars
        real(p2) :: num, denom
        !
        ! integer :: i, j

        ! Omega = zero
        ! do i = 2,4 ! gradQ, Q = [p u v w T]'
        !     do j = 1,3 ! if we really wanted to optimize this we could unroll this loop and remove diagonal terms...
        !         Wij   = half * (gradQ(j,i) - gradQ(i-1,j+1))
        !         Omega = Omega + two * ( Wij**2 )
        !     end do
        ! end do
        ! Omega = sqrt(Omega)
        uy = gradQ(iy, iu)
        uz = gradQ(iz, iu)
        vx = gradQ(ix, iv)
        vz = gradQ(iz, iv)
        wx = gradQ(ix, iw)
        wy = gradQ(iy, iw)
        
        ! The line below is equivalent to Omega
        Omega = sqrt( (uy - vx)**2 + (uz - wx)**2 + (vz - wy)**2 )
        ! dOmega = zero

        ! Alot of these functions have been optimized (note I never said I did a good job...)
        dChi = one / kvisc
        Chi = nut * dChi ! Chi = nut / kvisc

        num     = Chi**3
        denom   = (num + cv13) ! Chi**3 + cv1**3
        fv1     = num / denom    ! (Chi**3) / (Chi**3 + cv1**3)
        dfv1    = ( (three * dChi * Chi**2) * ( denom - num) ) / denom**2 ! Quotient rule: f' = g', f + C = g ==> f' * (g-f) / g**2
        
        denom   = (one + Chi * fv1)
        fv2     = one - Chi / denom
        dfv2    = - ( dChi - Chi * Chi * dfv1 ) / denom**2

        Sbar    = nut * fv2 / (KAPPA * distance)**2
        dSbar   = ( fv2 + zero * nut * dfv2) / (KAPPA * distance)**2 

        ! Equation (12) from ICCFD7-1902.  This is a modification of the standard implimentation that provides better numerical 
        ! stability. (See Note 1 from the NASA TMR for other potential treatments)
        if (Sbar >= -c2 * Omega) then
            Shat  = Omega + Sbar
            dShat = dSbar 
        else
            num   = c22 * Omega + c3 * Sbar
            denom = c3m2c2 * Omega - Sbar
            Shat  = Omega + Omega * num / denom
            dShat = dSbar * (c3 * Omega + num/denom) / denom
        endif
        if (Shat <= 1.e-010_p2) then
            shat = 1.e-010_p2
            dShat = zero
        endif

        ft2  = ct3 * exp(-ct4 * Chi**2)
        dft2 = ft2 * (-ct4 * 2 * Chi * dChi)

        denom  = Shat * (KAPPA * distance)**2
        r      = min(nut / denom, 10.0_p2)
        if (r == 10.0_p2) then
            dr = zero
        else
            dr = (denom - nut * dShat * (KAPPA * distance)**2) / denom**2
        endif

        g  = r + cw2 * (r**6 - r)
        dg = dr * ( one + cw2 * (six * r**5 - one) )

        num   = one + cw36
        denom = g**6 + cw36
        fw  = g * ( num / denom )**sixth ! float point exp. gross!
        dfw = dg * ( num / denom )**sixth - &
              g * sixth * ( six * num * g**5 ) / denom**2 / &
              ( num / denom )**fivesixth

        p     = cb1 * (one - ft2) * Shat
        prod  = p * nut
        pprm  = cb1 * ( (one-ft2)*dShat - dft2 * Shat)
        ! dprod = cb1 * ( ((one - dft2) * Shat * nut ) + (one-ft2) * (dShat * nut + Shat) )

        d     = ( cw1 * fw - (cb1 / KAPPA**2) * ft2 ) / distance**2
        ! dprm  = ( cw1 * dfw - cb1 * dft2 / KAPPA**2) * nut / distance**2 + d
        dprm  = d
        d     = d * nut
        dest  = d * nut
        ! ddest = ( ( cw1 * dfw - cb1 * dft2 / KAPPA**2) * (nut / distance)**2 ) + &
        !         ( cw1 * fw - (cb1 / KAPPA**2) * ft2 ) * two * nut / distance**2

        s1 = cb2 * dot_product(gradnut,gradnut) * iSIGMA 
        ! ds1 = zero

        ! if (prod > 1.0_p2) then
        !     write(*,*)
        ! endif
        source = prod - dest! + s1 
        
        ! Jacobian is treated using EQ (40) from https://doi.org/10.2514/6.1992-439 (multiplied by -1)
        ! jac = Pbar - Dbar = neg(prod - dest) + neg(dprod - ddest)*nut
        ! neg(x) = min(0,x)
        ! This treatment ensures diagonal dominance which improves convergence and stability of the linear solver.
        dsource = min(p - d, zero) + min(pprm - dprm, zero) * nut
        
        ! Trying a simpler approach
        ! dsource = min(dprod, zero) - max(ddest, zero)
        
    end subroutine sa_source


    
    ! pure function sa_prod(nut,ft2,strain_rate) result(source)

    !     use common , only : one

    !     real(p2), intent(in)  :: nut, ft2, strain_rate
    !     real(p2)              :: source

    !     source = cb1 * (one - ft2) * strain_rate * nut

    ! end function sa_prod

    ! pure function sa_dest(nut,distance,fw,ft2) result(dest)

    !     real(p2), intent(in) :: nut, distance, fw, ft2
    !     real(p2)             :: dest

    !     dest = ( cw1 * fw - (cb1 / KAPPA**2) * ft2 ) * (nut / distance)**2

    ! end function sa_dest


end module res_sa