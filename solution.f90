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
    real(p2), dimension(:)    , pointer :: ur2      ! U_R^2 for Weiss-Smith Low Mach Preconditioning
    

    real(p2), dimension(:)    , pointer :: dtau  !pseudo time step
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


    ! Jacobian type has to be placed here to avoid circular dependencies.
    type jacobian_type
        real(p2), dimension(5,5)                :: diag     ! diagonal blocks of Jacobian matrix
        real(p2), dimension(:,:,:), allocatable :: off_diag ! off-diagonal blocks
        real(p2), dimension(5,5)                :: diag_inv ! inverse of diagonal blocks
        real(p2), dimension(5)                  :: RHS      ! Right hand side (b) of the linear system
    end type jacobian_type

    public :: kth_nghbr_of_1, kth_nghbr_of_2
    integer, dimension(:), allocatable :: kth_nghbr_of_1
    integer, dimension(:), allocatable :: kth_nghbr_of_2

    public :: jac
    type(jacobian_type), dimension(:), allocatable :: jac ! jacobian array
    
    ! Force values
    real(p2) :: force_normalization

    contains

    subroutine allocate_solution_vars

        use common , only : p2, pi

        use grid , only : ncells, nnodes

        use config , only : accuracy_order, grad_method, lsq_stencil, solver_type, method_inv_flux, method_inv_jac, &
                            sideslip, aoa, lift, drag, turbulence_type

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

        if ( accuracy_order > 1 .OR. trim(turbulence_type) == 'laminar' ) then
            allocate( ccgradq(ndim,nq,ncells) )
            if (trim(grad_method) == 'lsq' .and. trim(lsq_stencil) == 'w_vertex') then
                allocate(  vgradq(ndim,nq,nnodes) )
            endif
        endif

        allocate( res(nq,ncells) )
        res = zero

        if ( trim(solver_type) == 'implicit') allocate(solution_update(nq,ncells))

        if ( trim(method_inv_flux) == 'roe_lm_w' ) then
            allocate(ur2(ncells))
        elseif (trim(method_inv_jac) == 'roe_lm_w' ) then
            write(*,*) "Weiss-Smith Low Mach Roe FDS can only be used for jacobians if it is also used for invicid flux."
            write(*,*) "method_inv_flux = ", trim(method_inv_flux)
            write(*,*) "method_inv_jac  = ", trim(method_inv_flux)
            write(*,*) "STOP"
            stop
        endif

        if ( lift ) then 
            ! These might even be correct :)
            vector_lift(1) = -sin(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
            vector_lift(2) = -sin(aoa*pi/180_p2)*sin(sideslip*pi/180_p2)
            vector_lift(3) =  cos(aoa*pi/180_p2)
        endif
        if ( drag ) then
            ! ditto:)
            vector_drag(1) =  cos(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
            vector_drag(2) =  cos(aoa*pi/180_p2)*sin(sideslip*pi/180_p2)
            vector_drag(3) =  sin(aoa*pi/180_p2)
        endif

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

    subroutine compute_uR2

        use common , only : p2, third, half, one, three_half

        use config , only : eps_weiss_smith, accuracy_order

        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp

        use grid   , only : ncells, cell

        use grid_statists , only : grid_spacing

        implicit none

        integer  :: icell
        real(p2) :: clength, dp, rho
        real(p2) :: C0, mu, T
        C0= sutherland_constant/reference_temp

        do icell = 1,ncells
            ! if (accuracy_order == 2) then !dp term
            !     ! clength = cell(icell)%vol**third
            !     dp = abs (  sqrt(ccgradq(1,1,i)**2 + ccgradq(2,1,i)**2 + ccgradq(3,1,i)**2) )
            !     rho = q(1,i) * gamma / q(5,i)
            ! endif
            ! ur2(i) = ( min( max( 0.001_p2, 0.1_p2 * sqrt(dp/rho),sqrt(q(2,i)**2 + q(3,i)**2 + q(4,i)**2) ), one) )**2
            clength = grid_spacing(icell)
            ! dp ~= dp/dx * dx ( just don't show this to a mathematician)
            dp = clength * abs (  sqrt(ccgradq(1,1,icell)**2 + ccgradq(2,1,icell)**2 + ccgradq(3,1,icell)**2) )
            T = q(5,icell)
            rho = q(1,icell) * gamma / T

            mu =  M_inf/Re_inf * (one + C0/T_inf) / (T + C0/T_inf)*T**(three_half)

            ur2(icell) = ( min( max( eps_weiss_smith * sqrt(dp/rho), &
                                     mu/rho/clength , &
                                     sqrt(q(2,icell)**2 + q(3,icell)**2 + q(4,icell)**2) ), &
                                     T ) )**2
        end do

    end subroutine compute_uR2

end module solution

module initialize

    implicit none

    public :: set_initial_solution
    public :: init_jacobian

    contains

    subroutine set_initial_solution

        use common , only : p2, one, pi, two

        use grid   , only : ncells

        use config , only : M_inf, aoa, sideslip, perturb_initial, random_perturb, solver_type, lift, drag, area_reference, &
                            high_ar_correction, method_inv_flux, method_inv_jac

        use solution

        use grid_statists , only : init_ar_array, compute_aspect_ratio, compute_grid_spacing

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
        
        if (trim(solver_type) == 'implicit' ) call init_jacobian
        
        force_normalization = two / ( rho_inf * area_reference *  M_inf**2 )

        if (high_ar_correction) then
            call init_ar_array
            call compute_aspect_ratio
        endif

        if (trim(method_inv_flux) == 'roe_lm_w' .OR. trim(method_inv_jac) == 'roe_lm_w' ) then
            call compute_grid_spacing
        endif
        
    end subroutine set_initial_solution

    
    
    subroutine init_jacobian

        use common          , only : p2

        use grid            , only : nfaces, face, cell, ncells

        use solution        , only : nq, jacobian_type, kth_nghbr_of_1, kth_nghbr_of_2, jac

        implicit none

        integer :: i, k
        integer :: c1, c2

        ! Create kth_nghbr arrays
        allocate(kth_nghbr_of_1(nfaces))
        allocate(kth_nghbr_of_2(nfaces))
        allocate(jac           (ncells))

        ! Define kth neighbor arrays
        face_nghbr_loop : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            ! loop over c1 neighbors to find c2
            do k = 1,cell(c1)%nnghbrs
                if ( c2 == cell(c1)%nghbr(k)) then
                    kth_nghbr_of_1(i) = k ! c2 is the kth neighbor of c1
                end if
            end do
            ! repeat for cell 2
            do k = 1,cell(c2)%nnghbrs
                if ( c1 == cell(c2)%nghbr(k)) then
                    kth_nghbr_of_2(i) = k
                end if
            end do
        end do face_nghbr_loop

        ! allocate jacobian off diagonal arrays
        do i = 1,ncells
            allocate(  jac(i)%off_diag(nq,nq,cell(i)%nnghbrs))
        end do

    end subroutine init_jacobian

end module initialize