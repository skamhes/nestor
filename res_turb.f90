module res_turb

    use common , only : p2
    
    implicit none

    public :: compute_residual_turb

    private
    
    contains

    subroutine compute_residual_turb

        use utils , only : iturb_type, TURB_SA

        use res_sa , only: compute_res_sa

        select case(iturb_type)
        case(TURB_SA)
            call compute_res_sa
        case default
            write(*,*) " Unsupported turbulence model. Stop"
            write(*,*) " res_turb.f90"
            stop
        end select
        

    end subroutine compute_residual_turb

end module res_turb