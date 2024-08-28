module ad_viscous_flux
    
    public :: visc_flux_internal_ddt
    public :: visc_flux_boundary_ddt

    private
    contains

    subroutine visc_flux_internal_ddt(q1,q2,gradq1,gradq2,n12,xc1,yc1,zc1,xc2,yc2,zc2,dFndQL,dFndQR)

        ! Face gradient terms computed using EQ. 14 in https://doi.org/10.2514/2.689 

        use common                  , only : p2, half, one, zero, three_half, two_third, four_third

        use solution                , only : gammamo, nq, ndim, T_inf ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use ad_operators

        implicit none

        real(p2), dimension(nq),      intent(in)  :: q1, q2
        real(p2), dimension(ndim,nq), intent(in)  :: gradq1, gradq2
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
                gradq_face(:,ivar) = gradq_face(:,ivar) + ( (qR_ddt(ivar) - qL_ddt(ivar)) - ddt_dot_product(gradq_face(:,ivar),ds,nq)) * dsds2
            end do

            ! This subroutine only handles computing the interface gradient.
            ! Once we have it we call the internal function
            call compute_visc_num_flux_ddt(qL_ddt,qR_ddt,gradq_face,n12,dFndQ)
            if (icell==1) then
                dFndQL = dFndQ
            else
                dFndQR = dFndQ
            endif
        end do jac_L_R
    end subroutine visc_flux_internal_ddt

    subroutine visc_flux_boundary_ddt(q1,q2,interface_grad_dummy,n12,dFndQL,dFndQR)

        use common                  , only : p2, half, one, zero, three_half, two_third, four_third

        use solution                , only : gammamo, nq, ndim
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use ad_operators

        implicit none

        real(p2), dimension(nq),      intent(in)  :: q1, q2
        real(p2), dimension(ndim,nq), intent(in)  :: interface_grad_dummy
        real(p2), dimension(ndim),    intent(in)  :: n12
        real(p2), dimension(nq,nq),   INTENT(OUT) :: dFndQL, dFndQR

        type(derivative_data_type_df5), dimension(ndim,nq)                   :: interface_grad
        type(derivative_data_type_df5), dimension(5) :: qL_ddt, qR_ddt
        real(p2),                       dimension(nq,nq)   :: dFndQ

        integer :: icell

        ! turn it into a ddt type (with zero df values).
        interface_grad = interface_grad_dummy

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

            ! This is just a wrapper function since we already have the interface gradient computed.
            call compute_visc_num_flux_ddt(qL_ddt,qR_ddt,interface_grad,n12,dFndQ)

            if (icell==1) then
                dFndQL = dFndQ
            else
                dFndQR = dFndQ
            endif
        end do jac_L_R
    end subroutine visc_flux_boundary_ddt

    subroutine compute_visc_num_flux_ddt(q1,q2,interface_grad,n12,dFdU)
        use common                  , only : p2, half, one, zero, three_half, two_third, four_third

        use solution                , only : gammamo, nq, ndim, T_inf ! w2u, nq
        
        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use ad_operators

        implicit none

        type(derivative_data_type_df5), dimension(nq),      intent(in)    :: q1, q2
        type(derivative_data_type_df5), dimension(ndim,nq), intent(in)    :: interface_grad
        real(p2), dimension(ndim),                          intent(in)    :: n12
        real(p2), dimension(nq,nq),                         intent(out)   :: dFdU

        ! Local Vars
        type(derivative_data_type_df5), dimension(nq) :: num_flux
        type(derivative_data_type_df5)                :: mu
        type(derivative_data_type_df5)                :: u, v, w, T
        type(derivative_data_type_df5)                :: C0
        type(derivative_data_type_df5)                :: tauxx, tauyy, tauzz !Viscous stresses: diagonal compontens
        type(derivative_data_type_df5)                :: tauxy, tauyz, tauzx !Viscous stresses: off-diagonal components
        type(derivative_data_type_df5)                :: tauyx, tauzy, tauxz !Viscous stresses: same as above by symmetry
        type(derivative_data_type_df5)                :: qx, qy, qz          !Heat flux components
        type(derivative_data_type_df5)                :: tauxn, tauyn, tauzn !Normal stresses
        type(derivative_data_type_df5)                :: qn                  !Normal heat flux
        type(derivative_data_type_df5), dimension(3)  :: grad_u, grad_v, grad_w   !Interface gradients of velocities
        type(derivative_data_type_df5), dimension(3)  :: grad_T

        integer :: i, j
        
        ! u = half * (q1(2)  + q1(2) ) ! u at the face
        ! v = half * (q1(3)  + q1(3) ) ! v at the face
        ! w = half * (q1(4)  + q1(4) ) ! w at the face
        T = half * (q1(nq) + q1(nq)) ! T at the face
        C0= sutherland_constant/reference_temp
        mu =  M_inf/Re_inf * (one + C0/T_inf) / (T + C0/T_inf)*T**(three_half)

        ! get_viscosity = scaling_factor * ( (one + ( C_0/Freestream_Temp ) )/(T + ( C_0/Freestream_Temp )) ) ** 1.5_p2
        if (ddt_isnan(mu)) then 
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

        do i = 1, 5
            do j = 1, 5
              !Flux derivative
               dFdU(i,j) = num_flux(i)%df(j)
            end do
        end do

    end subroutine compute_visc_num_flux_ddt

end module ad_viscous_flux