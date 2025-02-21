module turb_bc

    use common , only : p2

    implicit none
    
    contains

    subroutine sa_rhstate(nutL,bc_state_type, nutB)

        use common , only : p2, zero

        use turb,only : nut_inf

        use utils , only : BC_BACK_PRESSURE, BC_FARFIELD, BC_TANGENT, BC_VISC_STRONG

        implicit none

        ! Input 
        real(p2), intent(in) :: nutL
        integer , intent(in) :: bc_state_type

        ! Output
        real(p2), intent(out):: nutB

        select case(bc_state_type)
        case(BC_FARFIELD)
            nutB = nut_inf
        case(BC_TANGENT)
            nutB = nutL
        case(BC_VISC_STRONG) ! Adiabatic
            nutB = zero
        case(BC_BACK_PRESSURE)
            nutB = nut_inf ! This seems correct
        case default
            write(*,*) "Boundary condition=", bc_state_type ,"  not implemented."
            write(*,*) " --- Stop at get_right_state() in bc_states.f90..."
            stop
    end select
    end subroutine sa_rhstate

end module turb_bc