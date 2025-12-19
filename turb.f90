module turb

    use common , only : p2, three

    use utils  , only : iturb_model, iflow_type, FLOW_RANS, FLOW_DES, FLOW_LES, TURB_SA, &
                        isolver_type, SOLVER_IMPLICIT, SOLVER_GCR

    implicit none

    private

    ! VARS
    public turb_var
    public turb_jac
    public turb_res
    public ccgrad_turb_var, vgrad_turb_var
    public nturb
    public phi_turb
    public turb_jacobian_type
    public nut_inf
    public turb_res_norm, turb_res_norm_init
    public turb_update

    ! FUNCTIONS
    public allocate_rans
    public init_turb
    public calcmut

    real(p2), dimension(:,:)  , allocatable :: turb_var
    real(p2), dimension(:,:)  , allocatable :: turb_res
    real(p2), dimension(:)    , allocatable :: turb_update
    real(p2), dimension(:,:,:), allocatable :: ccgrad_turb_var, vgrad_turb_var
    integer                                 :: nturb
    real(p2), dimension(:)    , allocatable :: phi_turb
    real(p2), dimension(:)    , allocatable :: twsn, tdtau
    real(p2)                                :: nut_inf = three ! 3*mu/rho
    real(p2), dimension(7)                  :: turb_res_norm, turb_res_norm_init


    ! Jacobian type has to be placed here to avoid circular dependencies.
    type turb_jacobian_type
        real(p2)                                :: diag     ! diagonal blocks of Jacobian matrix
        real(p2), dimension(:), allocatable     :: off_diag ! off-diagonal blocks
        real(p2)                                :: diag_inv ! inverse of diagonal blocks
        real(p2)                                :: RHS      ! Right hand side (b) of the linear system
    end type turb_jacobian_type

    type(turb_jacobian_type), dimension(:,:), allocatable :: turb_jac

    contains

    subroutine allocate_rans
        
        use viscosity , only : compute_viscosity

        use grid      , only : ncells, nnodes, nfaces

        use solution_vars , only : kth_nghbr_of_1, kth_nghbr_of_2

        implicit none

        if (iflow_type > FLOW_RANS) then
            write(*,*) "DES and LES not yet supported"
            stop
        end if

        if (iturb_model == TURB_SA) then
            nturb = 1
        endif
        

        allocate(turb_var(ncells,nturb)) ! we're generally gonna be working through one variable at a time
        allocate(turb_res(ncells,nturb))
        allocate(turb_jac(ncells,nturb))
        
        allocate(turb_update(ncells))

        allocate(vgrad_turb_var( 3,nnodes,nturb))
        allocate(ccgrad_turb_var(3,ncells,nturb))

        allocate(phi_turb(ncells))

        allocate(twsn(ncells))
        allocate(tdtau(ncells))

        ! This isn't actually used except to prevent a runtime error.
        if(.not.allocated(kth_nghbr_of_1)) then
            allocate(kth_nghbr_of_1(nfaces))
            kth_nghbr_of_1 = 1
        endif
        if(.not.allocated(kth_nghbr_of_2)) then
            allocate(kth_nghbr_of_2(nfaces))
            kth_nghbr_of_2 = 1
        endif

    end subroutine allocate_rans

    subroutine init_turb

        use grid, only : ncells, cell

        use config    , only : turb_inf

        use solution_vars , only : mu_inf, rho_inf

        implicit none

        integer :: icell, it

        ! allocate jacobian off diagonal arrays
        do it = 1,nturb
            do icell = 1,ncells
                allocate( turb_jac(icell,it)%off_diag(cell(icell)%nnghbrs))
            end do
        end do

        ! Set freestream values
        if (iturb_model == TURB_SA) then
            nut_inf = turb_inf(nturb) * mu_inf / rho_inf
            turb_var(:,1) = nut_inf
        endif

    end subroutine init_turb

    function calcmut(q,mu,trbv) result(mut)

        use common , only : half

        use solution_vars , only : nq

        use utils , only : TURB_SA, iturb_model

        implicit none

        real(p2), dimension(nq),    intent(in) :: q
        real(p2),                   intent(in) :: mu
        real(p2), dimension(nturb), intent(in) :: trbv

        real(p2)                               :: mut

        select case(iturb_model)
        case(TURB_SA)
            mut = calcmut_SA(q,mu,trbv(1))
        case default
            write(*,*) "Error turb.f90 calcmut(). Unsupported turbulence model."
            stop
        end select  
        
    end function calcmut

    pure function calcmut_SA(q,mu,nut) result(mut)

    use common , only : half

    use solution_vars , only : nq, gamma, ip, iT

    use viscosity , only : compute_viscosity

    use sa_vars , only : cv13

    implicit none

    real(p2), dimension(nq), intent(in) :: q
    real(p2),                intent(in) :: mu
    real(p2),                intent(in) :: nut

    real(p2)                               :: mut

    ! Local Vars
    real(p2) :: rho
    real(p2) :: nu
    real(p2) :: chi3, fv1

    rho = q(ip)*gamma / q(iT)
    
    nu = mu / rho

    chi3 = ( nut / nu )**3

    fv1 = chi3 / ( chi3 + cv13 )

    mut = rho * nut * fv1
    
end function calcmut_SA

end module turb