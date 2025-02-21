module bc_states

    implicit none

    public

    contains

    subroutine get_right_state(qL,njk, bc_state_type, qcB)

        use common     , only : p2

        use utils , only : BC_BACK_PRESSURE, BC_FARFIELD, BC_TANGENT, BC_VISC_STRONG

        implicit none

        !Input
        real(p2), dimension(5),     intent(in) :: qL
        real(p2), dimension(3),     intent(in) :: njk
        integer ,           intent(in)         :: bc_state_type
        
        !output
        real(p2), dimension(5),    intent(out) :: qcB

        select case(bc_state_type)
            case(BC_FARFIELD)
                call freestream(qcB)
            case(BC_TANGENT)
                call slip_wall(qL,njk,qcB)
            case(BC_VISC_STRONG) ! Adiabatic
                call no_slip_wall(qL,qcB)
            case(BC_BACK_PRESSURE)
                call back_pressure(qL,qcB)
            case default
                write(*,*) "Boundary condition=", bc_state_type ,"  not implemented."
                write(*,*) " --- Stop at get_right_state() in bc_states.f90..."
                stop
        end select

    end subroutine get_right_state

    subroutine freestream(qb)
        use common      , only : p2
        use solution_vars    , only : p_inf, u_inf, v_inf, w_inf, T_inf
        implicit none

        real(p2),dimension(5), intent(out) :: qb

        qb(1) = p_inf
        qb(2) = u_inf
        qb(3) = v_inf
        qb(4) = w_inf
        qb(5) = T_inf

    end subroutine freestream

    subroutine back_pressure(qL,qcB)
        use common      , only : p2
        use solution_vars    , only : p_inf
        implicit none
        real(p2), dimension(5), intent( in) :: qL
        real(p2),dimension(5), intent(out) :: qcb

        qcb(1) = p_inf !<- Just fix the pressure.

        qcb(2:5) = qL(2:5)
        
    end subroutine back_pressure

    subroutine slip_wall(qL,njk,qcB)
        use common      , only : p2
        implicit none

        real(p2), dimension(5), intent( in) :: qL
        real(p2), dimension(3), intent( in) :: njk
        real(p2), dimension(5), intent(out) :: qcB

        real(p2) :: un
        
        un = qL(2)*njk(1) + qL(3)*njk(2) + qL(4)*njk(3)
        qcB = qL
        ! Ensure zero normal velocity on average:
        qcB(2) = qL(2) - un*njk(1)
        qcB(3) = qL(3) - un*njk(2)
        qcB(4) = qL(4) - un*njk(3)

    end subroutine slip_wall

    subroutine no_slip_wall(qL,qcB)
        ! no slip wall with zero heat flux (adiabatic condition)
        use common      ,   only : p2
        implicit none

        real(p2), dimension(5), intent( in) :: qL
        real(p2), dimension(5), intent(out) :: qcB
        
        ! un = wL(2)*njk(1) + wL(3)*njk(2) + wL(4)*njk(3)
        
        qcB(1) = qL(1)
        qcB(2:4) = -qL(2:4) ! half * (qL + qcB) = zero
        qcB(5) = qL(5)

    end subroutine no_slip_wall
end module bc_states