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
    public ccgrad_turb_var, vgrad_turb_var
    public nturb
    public phi_turb

    real(p2), dimension(:,:)  , allocatable :: turb_var
    real(p2), dimension(:,:)  , allocatable :: turb_res
    real(p2), dimension(:,:,:), allocatable :: ccgrad_turb_var, vgrad_turb_var
    integer                                 :: nturb
    real(p2), dimension(:)    , allocatable :: phi_turb

    ! Jacobian type has to be placed here to avoid circular dependencies.
    type turb_jacobian_type
        real(p2)                                :: diag     ! diagonal blocks of Jacobian matrix
        real(p2), dimension(:), allocatable     :: off_diag ! off-diagonal blocks
        real(p2)                                :: diag_inv ! inverse of diagonal blocks
        real(p2)                                :: RHS      ! Right hand side (b) of the linear system
    end type turb_jacobian_type

    type(turb_jacobian_type), dimension(:,:), allocatable :: turb_jac

    contains

    subroutine init_rans
        
        use viscosity , only : compute_viscosity

        use solution  , only : rho_inf, nut_inf, mu_inf

        use grid      , only : ncells, nnodes

        use config    , only : turb_inf

        use common    , only : zero

        implicit none

        if (iturb_type > TURB_RANS) then
            write(*,*) "DES and LES not yet supported"
            stop
        end if

        if (iturb_model == TURB_SA) then
            nturb = 1
            nut_inf = turb_inf(nturb) * mu_inf / rho_inf
        endif

        allocate(turb_var(ncells,nturb)) ! we're generally gonna be working through one variable at a time
        allocate(turb_res(ncells,nturb))
        allocate(turb_jac(ncells,nturb))

        allocate(vgrad_turb_var( 3,nnodes,nturb))
        allocate(ccgrad_turb_var(3,ncells,nturb))

        allocate(phi_turb(ncells))

    end subroutine init_rans

    subroutine init_turb_jacobian

        use grid, only : ncells, cell

        implicit none

        integer :: icell, it

        ! allocate jacobian off diagonal arrays
        do it = 1,nturb
            do icell = 1,ncells
                allocate( turb_jac(icell,it)%off_diag(cell(icell)%nnghbrs))
            end do
        end do

    end subroutine init_turb_jacobian

end module turb