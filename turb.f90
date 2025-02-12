module turb

    use common , only : p2

    use utils  , only : iturb_model, iturb_type, TURB_RANS, TURB_DES, TURB_LES, TURB_SA, &
                        isolver_type, SOLVER_IMPLICIT, SOLVER_GCR

    use solution,only: jacobian_type

    implicit none

    private

    public turb_var
    public turb_jac
    public turb_res
    public nturb
    public turb_inf

    real(p2), dimension(:,:)  , allocatable :: turb_var
    real(p2), dimension(:,:)  , allocatable :: turb_res
    integer                                 :: nturb

    type(jacobian_type), dimension(:,:), allocatable :: turb_jac

    contains

    subroutine init_rans
        
        use viscosity , only : compute_viscosity

        use solution  , only : T_inf, rho_inf

        use grid      , only : ncells

        use config    , only : turb_inf

        use common    , only : zero

        implicit none

        real(p2) :: mu_inf

        if (iturb_type > TURB_RANS) then
            write(*,*) "DES and LES not yet supported"
            stop
        end if

        if (iturb_model == TURB_SA) then
            nturb = 1
            if (turb_inf(nturb) > zero) then
                mu_inf = compute_viscosity(T_inf)
                turb_inf(nturb) = mu_inf / rho_inf
            endif
        endif

        allocate(turb_var(ncells,nturb)) ! we're generally gonna be working through one variable at a time
        allocate(turb_res(ncells,nturb))
        allocate(turb_jac(ncells,nturb))


        

    end subroutine init_rans

end module turb