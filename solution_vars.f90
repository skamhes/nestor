module solution_vars

    use common , only : p2, zero, one, three

    implicit none

    public

    !------------------------------------------
    !>> Solution data
    !------------------------------------------

    integer                             :: nq       ! # of eqns/solns (5 for 3D Euler/NS).
    integer                             :: ndim     ! # Number of dimesnions (3 for 3D (duh...))
    real(p2), dimension(:,:)  , pointer :: u        ! conservative variables at cells center
    real(p2), dimension(:,:)  , pointer :: w        ! primitive (rho, u, v, w, p) variables at cells center.
    real(p2), dimension(:,:)  , pointer :: q        ! alternate primitive variables (p, u, v, w, T) at cells centers.
    real(p2), dimension(:,:,:), pointer :: ccgradw  ! gradients of w at cell center.
    real(p2), dimension(:,:,:), pointer :: vgradw   ! gradients of w at vertices (nodes).
    real(p2), dimension(:,:,:), pointer :: ccgradq  ! gradients of q at cell center.
    real(p2), dimension(:,:,:), pointer :: vgradq   ! gradients of q at vertices (nodes).
    real(p2), dimension(:)    , pointer :: phi      ! limiter of ccgradq

    real(p2), dimension(:)    , pointer :: dtau  !pseudo time step
    real(p2)                            :: CFL_used ! used for GCR as CFL changes most iterations 
    real(p2), dimension(:)    , pointer :: wsn   !maximum eigenvalue at faces

    real(p2), dimension(:,:), pointer   :: res     !residual vector
    real(p2), dimension(5)              :: res_norm, res_norm_initial

    integer                             :: lrelax_sweeps_actual
    real(p2)                            :: lrelax_roc
    real(p2)                            :: roc

    real(p2), dimension(:,:), pointer   :: solution_update

    real(p2)                            :: force_drag
    real(p2)                            :: force_lift
    real(p2), dimension(3)              :: vector_drag
    real(p2), dimension(3)              :: vector_lift
    ! Note: I don't currently plan to use u and w (and gradw).  Instead my working 
    ! vars will be q.  But I'm keeping them as an option in case I change my mind 
    ! down the road.

    !------------------------------------------
    !>> Paramteres
    !------------------------------------------

    real(p2), parameter :: gamma = 1.4_p2 !Ratio of specific heats for air 
    ! For now this is a parameter (can't be changed).  In the future I may add support for different gamma values
    real(p2), parameter :: omgamma = one - gamma ! for convenience
    real(p2), parameter :: gammamo = gamma - one
    real(p2), parameter :: gmoinv  = one / (gammamo)

    !Free stream values: will be set in 'set set_initial_solution' for a given
    !free stream Mach number.
    real(p2) :: rho_inf = one
    real(p2) ::   u_inf = one
    real(p2) ::   v_inf = zero
    real(p2) ::   w_inf = zero
    real(p2) ::   p_inf = one/1.4_p2
    real(p2) ::   T_inf = one ! p_inf*gamma/rho_inf
    real(p2) ::   mu_inf= one ! will need to be overwritten

    ! Turb freestream values
    real(p2) ::   nut_inf = three ! 3*mu/rho

    !These data will be allocated for a given grid size, and filled in the
    !following subroutine: construct_ccfv_data.

    !------------------------------------------------------------------------------------
    ! End of the data used to implement a cell-centered FV method.
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------

    ! because I just had a bug where i used the wrong index.......
    integer, parameter :: ip = 1
    integer, parameter :: iu = 2
    integer, parameter :: iv = 3
    integer, parameter :: iw = 4
    integer, parameter :: iT = 5

end module solution_vars