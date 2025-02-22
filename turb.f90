module turb

    use common , only : p2, three

    use utils  , only : iturb_model, iturb_type, TURB_RANS, TURB_DES, TURB_LES, TURB_SA, &
                        isolver_type, SOLVER_IMPLICIT, SOLVER_GCR

    implicit none

    private

    public turb_var
    public turb_jac
    public turb_res
    public ccgrad_turb_var, vgrad_turb_var
    public nturb
    public phi_turb
    public turb_jacobian_type
    public nut_inf
    public allocate_rans
    public init_turb
    public turb_res_norm
    public turb_update

    real(p2), dimension(:,:)  , allocatable :: turb_var
    real(p2), dimension(:,:)  , allocatable :: turb_res
    real(p2), dimension(:)    , allocatable :: turb_update
    real(p2), dimension(:,:,:), allocatable :: ccgrad_turb_var, vgrad_turb_var
    integer                                 :: nturb
    real(p2), dimension(:)    , allocatable :: phi_turb
    real(p2)                                :: nut_inf = three ! 3*mu/rho
    real(p2), dimension(7)                  :: turb_res_norm

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

        use grid      , only : ncells, nnodes

        implicit none

        if (iturb_type > TURB_RANS) then
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

end module turb