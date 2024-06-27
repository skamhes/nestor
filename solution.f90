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
    real(p2), dimension(:,:)  , pointer :: u        ! conservative variables at cells center
    real(p2), dimension(:,:)  , pointer :: w        ! primitive (rho, u, v, w, p) variables at cells center.
    real(p2), dimension(:,:)  , pointer :: q        ! alternate primitive variables (p, u, v, w, T) at cells centers.
    real(p2), dimension(:,:,:), pointer :: ccgradw  ! gradients of w at cell center.
    real(p2), dimension(:,:,:), pointer :: vgradw   ! gradients of w at vertices (nodes).
    real(p2), dimension(:,:,:), pointer :: ccgradq  ! gradients of q at cell center.
    real(p2), dimension(:,:,:), pointer :: vgradq   ! gradients of q at vertices (nodes).

    real(p2), dimension(:)    , pointer :: dtau  !pseudo time step
    real(p2), dimension(:)    , pointer :: wsn   !maximum eigenvalue at faces

    real(p2), dimension(:,:), pointer :: res     !residual vector
    real(p2), dimension(5)            :: res_norm, res_norm_initial

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


    !These data will be allocated for a given grid size, and filled in the
    !following subroutine: construct_ccfv_data.

    !------------------------------------------------------------------------------------
    ! End of the data used to implement a cell-centered FV method.
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------

    contains

    subroutine allocate_solution_vars

        use grid , only : ncells, nnodes

        use config , only : accuracy_order, grad_method, lsq_stencil

        implicit none

        nq = 5

        ! initialize
        allocate( q(nq,ncells) )
        q = zero

        allocate( dtau(ncells) )
        allocate( wsn( ncells) )
        dtau = zero
        wsn = zero

        if ( accuracy_order > 1 ) then
            allocate( ccgradq(3,nq,ncells) )
            if (trim(grad_method) == 'lsq' .and. trim(lsq_stencil) == 'wvertex') then
                allocate(  vgradq(3,nq,nnodes) )
            endif
        endif

        allocate( res(nq,ncells) )
        res = zero

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
        q_out(5) = q_out(1)*gamma/w_in(1)
  
    end function w2q

end module solution

module initialize

    implicit none

    public set_initial_solution

    contains

    subroutine set_initial_solution

        use common , only : p2, one, pi

        use grid   , only : ncells

        use config , only : M_inf, aoa, sideslip, perturb_initial

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

        cell_loop : do i = 1,ncells
            q(:,i) = q_init
        end do cell_loop
        
    end subroutine set_initial_solution
end module initialize