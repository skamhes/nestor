module ad_viscous_flux
    
    public :: visc_flux_internal_ddt
    public :: visc_flux_boundary_ddt

    private
    contains

    subroutine visc_flux_internal_ddt(q1,q2,gradq1,gradq2,trb1,trb2,n12,xc1,yc1,zc1,xc2,yc2,zc2,dFndQL,dFndQR)

        ! Face gradient terms computed using EQ. 14 in https://doi.org/10.2514/2.689 

        use common                  , only : p2, half

        use solution_vars           , only : nq, ndim ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use ad_operators

        use turb                    , only : nturb

        implicit none

        real(p2), dimension(nq),      intent(in)  :: q1, q2
        real(p2), dimension(ndim,nq), intent(in)  :: gradq1, gradq2
        real(p2), dimension(nturb),   intent(in)  :: trb1, trb2
        real(p2), dimension(ndim),    intent(in)  :: n12               ! Unit area vector (from c1 to c2)
        real(p2),                     intent(in)  :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                     intent(in)  :: xc2, yc2, zc2     ! Right cell centroid
        real(p2), dimension(nq,nq),   INTENT(OUT) :: dFndQL, dFndQR

        ! Local Vars
        type(derivative_data_type_df5), dimension(5) :: qL_ddt, qR_ddt
        type(derivative_data_type_df5), dimension(ndim,nq) :: gradq_face
        real(p2),                       dimension(ndim)    :: ds,  dsds2
        real(p2),                       dimension(nq,nq)   :: dFndQ
        
        integer                      :: icell, ivar

        jac_L_R : do icell = 1,2
            qL_ddt = q1
            qR_ddt = q2
            if (icell == 1) then
                ! Using derivative for uL_ddt
                call ddt_seed(qL_ddt)
            else ! icell = 2
                ! Using derivative for uR_ddt
                call ddt_seed(qR_ddt)
            end if
            ! Calculate the face gradients
            ds = (/xc2-xc1, yc2-yc1, zc2-zc1/) ! vector pointing from center of cell 1 to cell 2
            dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds(:)/ds**2

            ! Equation 14
            do ivar = 1,nq
                gradq_face(:,ivar) = half * (gradq1(:,ivar) + gradq2(:,ivar))
                gradq_face(:,ivar) = gradq_face(:,ivar) + & 
                                     ( (qR_ddt(ivar) - qL_ddt(ivar)) - ddt_dot_product(gradq_face(:,ivar),ds,ndim)) * dsds2
            end do

            ! This subroutine only handles computing the interface gradient.
            ! Once we have it we call the internal function
            call compute_visc_num_flux_ddt(qL_ddt,qR_ddt,trb1,trb2,gradq_face,n12,dFndQ)
            if (icell==1) then
                dFndQL = dFndQ
            else
                dFndQR = dFndQ
            endif
        end do jac_L_R
    end subroutine visc_flux_internal_ddt

    subroutine visc_flux_boundary_ddt(q1,qb,interface_grad_dummy,trb1,trb2,n12,xc1,yc1,zc1,xf2,yf2,zf2,dFndQL,dFndQR)

        use common                  , only : p2, half

        use solution_vars           , only : nq, ndim
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use ad_operators

        use turb                    , only : nturb

        implicit none

        real(p2), dimension(nq),      intent(in)  :: q1, qb
        real(p2), dimension(nturb),   intent(in)  :: trb1, trb2
        real(p2), dimension(ndim,nq), intent(in)  :: interface_grad_dummy
        real(p2), dimension(ndim),    intent(in)  :: n12
        real(p2),                     intent(in)  :: xc1, yc1, zc1     ! Left cell centroid
        real(p2),                     intent(in)  :: xf2, yf2, zf2     ! Boundary face centroid
        real(p2), dimension(nq,nq),   INTENT(OUT) :: dFndQL, dFndQR

        ! Local Vars
        type(derivative_data_type_df5), dimension(ndim,nq) :: gradq_face
        type(derivative_data_type_df5), dimension(5)       :: qL_ddt, qR_ddt
        real(p2),                       dimension(nq,nq)   :: dFndQ
        real(p2), dimension(ndim)                          :: ds,  dsds2
        
        integer :: icell, ivar

        ds = (/xf2-xc1, yf2-yc1, zf2-zc1/) ! vector pointing from center of cell 1 to cell 2
        dsds2 = ds/(ds(1)**2 + ds(2)**2 + ds(3)**2) ! ds(:)/ds**2

        jac_L_R : do icell = 1,2
            qL_ddt = q1
            qR_ddt = qb
            if (icell == 1) then
                ! Using derivative for uL_ddt
                call ddt_seed(qL_ddt)
            else ! icell = 2
                ! Using derivative for uR_ddt
                call ddt_seed(qR_ddt)
            end if

            gradq_face = interface_grad_dummy

            ! Equation 14
            do ivar = 1,nq
                gradq_face(:,ivar) = gradq_face(:,ivar) + ( half * (qR_ddt(ivar) - qL_ddt(ivar)) &
                                     - ddt_dot_product(gradq_face(:,ivar),ds,ndim)) * dsds2
            end do

            ! This is just a wrapper function since we already have the interface gradient computed.
            call compute_visc_num_flux_ddt(qL_ddt,qR_ddt,trb1,trb2,gradq_face,n12,dFndQ)

            if (icell==1) then
                dFndQL = dFndQ
            else
                dFndQR = dFndQ
            endif
        end do jac_L_R
    end subroutine visc_flux_boundary_ddt

    subroutine compute_visc_num_flux_ddt(q1,q2,trb1,trb2,interface_grad,n12,dFdU)
        use common                  , only : p2, half, one, zero, three_half, two_third, four_third, ix, iy, iz

        use solution_vars           , only : gammamo, nq, ndim, iu, iv, iw, iT ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp, pr_t

        use ad_operators

        use viscosity               , only : compute_viscosity_ddt

        use turb                    , only : nturb, calcmut

        implicit none

        type(derivative_data_type_df5), dimension(nq),      intent(in)    :: q1, q2
        real(p2), dimension(nturb),   intent(in)                          :: trb1, trb2
        type(derivative_data_type_df5), dimension(ndim,nq), intent(in)    :: interface_grad
        real(p2), dimension(ndim),                          intent(in)    :: n12
        real(p2), dimension(nq,nq),                         intent(out)   :: dFdU

        ! Local Vars
        type(derivative_data_type_df5), dimension(nq) :: num_flux, qf
        type(derivative_data_type_df5)                :: mu, mu_effective
        type(derivative_data_type_df5)                :: p, u, v, w, T
        type(derivative_data_type_df5)                :: tauxx, tauyy, tauzz !Viscous stresses: diagonal compontens
        type(derivative_data_type_df5)                :: tauxy, tauyz, tauzx !Viscous stresses: off-diagonal components
        type(derivative_data_type_df5)                :: tauyx, tauzy, tauxz !Viscous stresses: same as above by symmetry
        type(derivative_data_type_df5)                :: qx, qy, qz          !Heat flux components
        type(derivative_data_type_df5)                :: tauxn, tauyn, tauzn !Normal stresses
        type(derivative_data_type_df5)                :: qn                  !Normal heat flux
        type(derivative_data_type_df5), dimension(3)  :: grad_u, grad_v, grad_w   !Interface gradients of velocities
        type(derivative_data_type_df5), dimension(3)  :: grad_T

        real(p2), dimension(nturb)                    :: trb
        real(p2)                                      :: mut ! not contributing to the Jacobian for now

        integer :: i, j
        
        qf = half * (q1 + q2)

        p = qf(1)  ! p at the face
        u = qf(2)  ! u at the face
        v = qf(3)  ! v at the face
        w = qf(4)  ! w at the face
        T = qf(5)  ! T at the face
        
        mu = compute_viscosity_ddt(T) ! for now we aren't counting this term towards the jacobian

        trb = half * (trb1 + trb2)

        mut = calcmut((/p%f,u%f,v%f,w%f,T%f/),mu%f,trb)

        mu_effective = mu + mut

#ifdef NANCHECK
        ! get_viscosity = scaling_factor * ( (one + ( C_0/Freestream_Temp ) )/(T + ( C_0/Freestream_Temp )) ) ** 1.5_p2
        if (ddt_isnan(mu_effective)) then 
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

        do i = 1, 5
            do j = 1, 5
              !Flux derivative
               dFdU(i,j) = num_flux(i)%df(j)
            end do
        end do

    end subroutine compute_visc_num_flux_ddt

end module ad_viscous_flux