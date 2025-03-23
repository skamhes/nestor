module turb_bc

    use common , only : p2

    implicit none
    
    contains

    subroutine turb_rhstate(turbvarL,bc_state_type, turbvarB)

        use utils , only : iturb_model, TURB_SA

        use turb , only : nturb

        implicit none

        ! Input 
        real(p2), dimension(:), intent(in) :: turbvarL
        integer ,               intent(in) :: bc_state_type

        ! Output
        real(p2), dimension(:), intent(out):: turbvarB

        select case(iturb_model)
        case(TURB_SA)
            call sa_rhstate(turbvarL(1),bc_state_type,turbvarB(1))
        case default
            write(*,*) " Unsupported turbulence model. Stop"
            write(*,*) " res_turb.f90"
            stop
        end select

    end subroutine turb_rhstate

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