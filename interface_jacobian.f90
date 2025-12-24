module interface_jacobian

    implicit none

    public :: interface_jac

    interface interface_jac
        module procedure interface_jac_std
        module procedure interface_jac_ws
    end interface interface_jac
    
    contains

    subroutine interface_jac_std(qj, qk, njk, dFndqL, dFndqR)

        use common              , only : p2, zero, one, half

        use ad_operators        ! all

        use utils               , only : imethod_inv_jac, IJAC_ROE, IJAC_HLL, IJAC_RHLL, IJAC_RUSANOV
        
        use ad_inviscid_flux    , only :      roe_ddt, &
                                            rusanov_ddt, &
                                                hll_ddt, &
                                               rhll_ddt
                                     
        use solution            , only : q2u

        implicit none

        real(p2), dimension(5), intent(in) :: qj, qk  ! w from cell(j) and neighbor (k)
        real(p2), dimension(3), intent(in) :: njk

        real(p2), dimension(5,5), intent(out) :: dFndqL, dFndqR

        ! Local vavrs
        real(p2), dimension(5,5)    :: dfndq
        real(p2), dimension(5)      :: dummy5
        real(p2)                    :: wsn
        
        integer :: i
        type(derivative_data_type_df5), dimension(5) :: uL_ddt, uR_ddt, qL_ddt, qR_ddt

        jac_L_R : do i = 1,2
            qL_ddt = qj
            qR_ddt = qk
            if (i == 1) then
                ! Using derivative for uL_ddt
                call ddt_seed(qL_ddt)
            else ! i = 2
                ! Using derivative for uR_ddt
                call ddt_seed(qR_ddt)
            end if
            select case(imethod_inv_jac)
            !------------------------------------------------------------
            !  (1) Roe flux
            !------------------------------------------------------------
            case(IJAC_ROE)
                ! call roe_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndq,wsn)        
                call roe_ddt(qL_ddt,qR_ddt,njk,dfndq)

            ! !------------------------------------------------------------
            ! !  (2) Rusanov flux
            ! !------------------------------------------------------------
            ! case(IJAC_RUSANOV)
            !     call rusanov_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndq,wsn)
            ! !------------------------------------------------------------
            ! !  (3) HLL flux
            ! !------------------------------------------------------------
            ! case(IJAC_HLL)
            !     call hll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndq,wsn)
            ! !------------------------------------------------------------
            ! !  (4) RHLL flux: the last argumant -> exact_jac = .false.
            ! !                  so that the jac = a1*HLL_jac+a2*Roe_jac
            ! !                   with a1 and a2 not differentiated.
            ! !------------------------------------------------------------
            ! case(IJAC_RHLL)
            !     call rhll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndq,wsn,.false.)
            !------------------------------------------------------------
            !  Others...
            !------------------------------------------------------------
            case default
                write(*,*) " Invalid input for inviscid_jac = ", imethod_inv_jac
                write(*,*) " Choose roe or rhll, and try again. interface_jacobian.f90"
                write(*,*) " ... Stop."
                stop
            end select

            if (i==1) then
                dFndqL = dfndq
            else
                dFndqR = dfndq
            endif
        end do jac_L_R

    end subroutine interface_jac_std

    subroutine interface_jac_ws(qj, qk, njk, ur2j, ur2k, dFnduL, dFnduR)

        use common              , only : p2, zero, one, half

        use ad_operators        ! all

        use config              , only : method_inv_jac

        use utils               , only : imethod_inv_jac, IJAC_ROE_LM
        
        use ad_inviscid_flux    , only :      roe_lm_w_ddt
                                     
        use solution            , only : q2u

        implicit none

        real(p2), dimension(5), intent(in) :: qj, qk  ! w from cell(j) and neighbor (k)
        real(p2), dimension(3), intent(in) :: njk
        real(p2),               intent(in) :: ur2j, ur2k  ! reference velocity from cell(j) and neighbor (k)

        real(p2), dimension(5,5), intent(out) :: dFnduL, dFnduR

        ! Local vavrs
        real(p2), dimension(5,5)    :: dfndu
        real(p2), dimension(5)      :: dummy5
        real(p2)                    :: wsn
        
        integer :: i
        type(derivative_data_type_df5), dimension(5) :: uL_ddt, uR_ddt, qL_ddt, qR_ddt

        jac_L_R : do i = 1,2
            qL_ddt = qj
            qR_ddt = qk
            if (i == 1) then
                ! Using derivative for uL_ddt
                call ddt_seed(qL_ddt)
                uL_ddt = q2u_ddt(qL_ddt)
                uR_ddt = q2u_ddt(qR_ddt)
            else ! i = 2
                ! Using derivative for uR_ddt
                call ddt_seed(qR_ddt)
                uL_ddt = q2u_ddt(qL_ddt)
                uR_ddt = q2u_ddt(qR_ddt)
            end if
            !------------------------------------------------------------
            !  (1) Roe flux
            !------------------------------------------------------------
            if(imethod_inv_jac == IJAC_ROE_LM) then
                    
                call roe_lm_w_ddt(uL_ddt,uR_ddt,njk,ur2j,ur2k,dummy5,dfndu,wsn)

            ! !------------------------------------------------------------
            ! !  (2) Rusanov flux
            ! !------------------------------------------------------------
            ! case(IJAC_RUSANOV)
            !     call rusanov_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndq,wsn)
            ! !------------------------------------------------------------
            ! !  (3) HLL flux
            ! !------------------------------------------------------------
            ! case(IJAC_HLL)
            !     call hll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndq,wsn)
            ! !------------------------------------------------------------
            ! !  (4) RHLL flux: the last argumant -> exact_jac = .false.
            ! !                  so that the jac = a1*HLL_jac+a2*Roe_jac
            ! !                   with a1 and a2 not differentiated.
            ! !------------------------------------------------------------
            ! case(IJAC_RHLL)
            !     call rhll_ddt(uL_ddt,uR_ddt,njk, dummy5,dfndq,wsn,.false.)
            !------------------------------------------------------------
            !  Others...
            !------------------------------------------------------------
            else
                write(*,*) " Invalid input for inviscid_jac = ", trim(method_inv_jac)
                write(*,*) " Choose roe or rhll, and try again."
                write(*,*) " ... Stop."
                stop
            endif
            if (i==1) then
                dFnduL = dfndu
            else
                dFnduR = dfndu
            endif
        end do jac_L_R

    end subroutine interface_jac_ws 

    !********************************************************************************
    ! Compute U from W (ddt version)
    !
    ! ------------------------------------------------------------------------------
    !  Input:  q =    primitive variables (  p,     u,     v,     w,     T)
    ! Output:  u = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
    ! ------------------------------------------------------------------------------
    !
    ! Note: rho*E = p/(gamma-1) + rho*0.5*(u^2 + v^2 + w^2)
    !       rho   = p/(gamma*T)
    !********************************************************************************
    function q2u_ddt(q_in) result(u_out)

        use common              , only : p2, one, half
        use solution            , only : gamma, gammamo, gmoinv
        use ad_operators

        implicit none
    
        type(derivative_data_type_df5), dimension(5), intent(in) :: q_in ! input
        type(derivative_data_type_df5), dimension(5)             :: u_out !output
    
        u_out(1) = q_in(1)*gamma / q_in(5)
        u_out(2) = u_out(1)*q_in(2)
        u_out(3) = u_out(1)*q_in(3)
        u_out(4) = u_out(1)*q_in(4)
        u_out(5) = q_in(1)*gmoinv+half*u_out(1)*(q_in(2)*q_in(2)+q_in(3)*q_in(3)+q_in(4)*q_in(4))
    
    end function q2u_ddt

end module interface_jacobian