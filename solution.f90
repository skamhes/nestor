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

    ! because I just had a bug where i used the wrong index.......
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
    
    ! ! Jacobian Free Newton-Krylov Variables
    ! real(p2), dimension(:,:,:), pointer :: gcr_precond_correction
    ! real(p2), dimension(:,:),   pointer :: gcr_final_correction
    ! real(p2), dimension(:,:,:), pointer :: gcr_search_direction
    real(p2) :: inv_ncells ! 1/ncells/nq
    real(p2) :: nl_reduction
    integer  :: n_projections                   




    ! Force values
    real(p2) :: force_normalization

    contains

    subroutine allocate_solution_vars

        use common , only : p2, pi

        use grid , only : ncells, nnodes

        use config , only : accuracy_order, grad_method, lsq_stencil, solver_type, lift, drag, aoa, sideslip, turbulence_type,&
                            gcr_max_projections

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

        if ( trim(solver_type) == 'implicit' .OR. trim(solver_type) == 'gcr') then 
            allocate(solution_update(nq,ncells))
        ! else if ( trim(solver_type) == 'gcr') then
        !     allocate(gcr_final_correction(nq,ncells))
        !     allocate(gcr_precond_correction(nq,ncells,gcr_max_projections))
        !     allocate(gcr_search_direction(  nq,ncells,gcr_max_projections))
        endif

        if ( lift ) then 
            ! These might even be correct :)
            vector_lift(1) = -sin(aoa*pi/180.0_p2)*cos(sideslip*pi/180.0_p2)
            vector_lift(2) = -sin(aoa*pi/180.0_p2)*sin(sideslip*pi/180.0_p2)
            vector_lift(3) =  cos(aoa*pi/180.0_p2)
        endif
        if ( drag ) then
            ! ditto:)
            vector_drag(1) =  cos(aoa*pi/180.0_p2)*cos(sideslip*pi/180.0_p2)
            vector_drag(2) =  cos(aoa*pi/180.0_p2)*sin(sideslip*pi/180.0_p2)
            vector_drag(3) =  sin(aoa*pi/180.0_p2)
        endif

        inv_ncells = one / real(ncells*nq,p2) 

    end subroutine allocate_solution_vars

    subroutine compute_local_time_step_dtau

        use common                  , only : half, p2
        use grid                    , only : ncells, cell
        use config                  , only : CFL, high_ar_correction
        use grid_statists           , only : cell_aspect_ratio

        implicit none

        integer :: i
        ! real(p2), dimension(ncells) :: viscous_dtau

        cell_loop : do i = 1,ncells
            dtau(i) = CFL * cell(i)%vol/( half * wsn(i) )
            if (high_ar_correction) dtau(i) = dtau(i) * cell_aspect_ratio(i)
        end do cell_loop
    end subroutine compute_local_time_step_dtau

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

    pure function compute_primative_jacobian(qi) result(preconditioner)

        ! This function computes the Jacobian DU/DQ where U is the vector of conserved variables [rho rhoU rhoV rhoW rhoE] and W is
        ! the vector of primitive variables [p U V W T].  This is also written in a way that can allow implementation of Weiss-Smith
        ! Preconditioning

        use common , only : p2, half

        implicit none

        real(p2), dimension(5),  intent(in) :: qi
        real(p2), dimension(5,5)            :: preconditioner
        
        real(p2), dimension(5,5) :: dwdq, pre_inv
        real(p2) :: H, rho_p, rho_T, theta, rho, uR2inv


        H = ((qi(iT))**2)*gmoinv + half * ( qi(iu)**2 + qi(iv)**2 + qi(iw)**2 )
        rho_p = gamma/qi(5)
        rho_T = - (qi(ip)*gamma)/(qi(iT)**2)
        rho = qi(ip)*gamma/qi(iT)
        UR2inv = one ! will be 1/uR2(i)
        theta = (UR2inv) - rho_T*(gammamo)/(rho)
        
        ! Note transposing this assignment would likely be marginally faster if slightly less easy to read
        preconditioner(1,:) = (/ theta,         zero,       zero,       zero,       rho_T                   /)
        preconditioner(2,:) = (/ theta*qi(iu),  rho,        zero,       zero,       rho_T*qi(iu)            /)
        preconditioner(3,:) = (/ theta*qi(iv),  zero,       rho,        zero,       rho_T*qi(iv)            /)
        preconditioner(4,:) = (/ theta*qi(iw),  zero,       zero,       rho,        rho_T*qi(iw)            /)
        preconditioner(5,:) = (/ theta*H-one,  rho*qi(iu), rho*qi(iv), rho*qi(iw), rho_T*H + rho/(gamma-one)/)

    end function compute_primative_jacobian

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
                            high_ar_correction

        use solution

        use grid_statists , only : init_ar_array, compute_aspect_ratio

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
        
        if (trim(solver_type) == 'implicit' .OR. trim(solver_type) == 'gcr' ) call init_jacobian
        
        force_normalization = two / ( rho_inf * area_reference *  M_inf**2 )

        if (high_ar_correction) then
            call init_ar_array
            call compute_aspect_ratio
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