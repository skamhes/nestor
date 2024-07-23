module solution

    ! This module contains all of the arrays and variables 
    ! needed by any of the solver modules
    use common , only : p2, one, zero
    
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
    real(p2), dimension(:)    , pointer :: wsn   !maximum eigenvalue at faces

    real(p2), dimension(:,:), pointer   :: res     !residual vector
    real(p2), dimension(5)              :: res_norm, res_norm_initial

    integer                             :: lrelax_sweeps_actual
    real(p2)                            :: lrelax_roc
    real(p2)                            :: roc

    real(p2), dimension(:,:), pointer   :: solution_update
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

    !These data will be allocated for a given grid size, and filled in the
    !following subroutine: construct_ccfv_data.

    !------------------------------------------------------------------------------------
    ! End of the data used to implement a cell-centered FV method.
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------

    ! because I just had a bug where i used the wronog index.......
    integer, parameter :: ip = 1
    integer, parameter :: iu = 2
    integer, parameter :: iv = 3
    integer, parameter :: iw = 4
    integer, parameter :: iT = 5

    contains

    subroutine allocate_solution_vars

        use grid , only : ncells, nnodes

        use config , only : accuracy_order, grad_method, lsq_stencil, solver_type

        implicit none

        nq = 5
        ndim = 3

        ! initialize
        allocate( q(nq,ncells) )
        q = zero

        allocate( dtau(ncells) )
        allocate( wsn( ncells) )
        dtau = zero
        wsn = zero

        if ( accuracy_order > 1 ) then
            allocate( ccgradq(ndim,nq,ncells) )
            if (trim(grad_method) == 'lsq' .and. trim(lsq_stencil) == 'w_vertex') then
                allocate(  vgradq(ndim,nq,nnodes) )
            endif
        endif

        allocate( res(nq,ncells) )
        res = zero

        if ( trim(solver_type) == 'implicit') allocate(solution_update(nq,ncells))

    end subroutine allocate_solution_vars

    !********************************************************************************
    ! Compute Q from W
    !
    ! ------------------------------------------------------------------------------
    !  Input:  u = conservative variables (rho, u, v, w, p)
    ! Output:  q =    primitive variables (  p, u, v, w, T)
    ! ------------------------------------------------------------------------------
    !
    ! Note:    E = p/(gamma-1)/rho + 0.5*(u^2+v^2+w^2)
    !       -> p = (gamma-1)*rho*E-0.5*rho*(u^2+v^2w^2)
    !       -> T = p*gamma/rho     
    ! 
    !********************************************************************************
    function w2q(w_in) result(q_out)
  
    use common, only : p2
    
    implicit none
  
    real(p2), dimension(5), intent(in) :: w_in ! input
    real(p2), dimension(5)             :: q_out !output
  
        q_out(2) = w_in(2)
        q_out(3) = w_in(3)
        q_out(4) = w_in(4)
        
        q_out(1) = w_in(5)
        ! T      = p       *gamma/rho
        q_out(5) = q_out(1)*gamma/w_in(1)
  
    end function w2q

    !********************************************************************************
    ! Compute Q from U
    !
    ! ------------------------------------------------------------------------------
    !  Input:  Q =    primitive variables (  p,     u,     v,     w,     T)
    ! Output:  U = conservative variables (rho, rho*u, rho*v, rho*w, rho*E)
    ! ------------------------------------------------------------------------------
    !
    ! Note: rho*E = p/(gamma-1) + rho*0.5*(u^2 + v^2 + w^2)
    !       rho   = p*gamma/T
    !********************************************************************************
    function q2u(q_in) result(u_out)

        use common, only : p2, half
        
        implicit none
    
        real(p2), dimension(5), intent(in) :: q_in ! input
        real(p2), dimension(5)             :: u_out !output
        
        ! rho    = p      *gamma / T
        u_out(1) = q_in(1)*gamma / q_in(5)
        u_out(2) = u_out(1)*q_in(2)
        u_out(3) = u_out(1)*q_in(3)
        u_out(4) = u_out(1)*q_in(4)
        u_out(5) = q_in(1)*gmoinv + half*u_out(1)*(q_in(iu)**2 + q_in(iv)**2 + q_in(iw)**2)
    
    end function q2u

end module solution

module initialize

    implicit none

    public set_initial_solution

    contains

    subroutine set_initial_solution

        use common , only : p2, one, pi

        use grid   , only : ncells

        use config , only : M_inf, aoa, sideslip, perturb_initial, random_perturb

        use solution

        implicit none

        integer                 :: i
        real(p2), dimension(5)  :: q_init

        ! Set the free stream values
        rho_inf = one
        u_inf = M_inf*cos(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        v_inf = M_inf*sin(sideslip*pi/180)
        w_inf = M_inf*sin(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        p_inf = one/gamma

        q_init = w2q( (/rho_inf,u_inf,v_inf,w_inf,p_inf/) )
        if ( perturb_initial )  then 
            q_init(2:4) = (/ 0.2_p2, 0.1_p2, 0.15_p2 /)
        end if

        cell_loop : do i = 1,ncells
        q(:,i) = q_init
            if ( perturb_initial .and. random_perturb )  then 
                q(2:4,i) = q(2:4,i) * rand(0)
            endif
        end do cell_loop
        
    end subroutine set_initial_solution
end module initialize