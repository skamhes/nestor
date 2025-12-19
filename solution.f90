module solution

    contains

    subroutine define_problem

        ! This just avoids some circular dependencies for now.  But later on we can use it to define additional solvers that we may 
        ! need.

        nq = 5
        ndim = 3

    end subroutine define_problem

    subroutine allocate_solution_vars

        use common , only : p2, pi

        use grid , only : ncells, nnodes

        use config , only : accuracy_order, grad_method, lsq_stencil, lift, drag, aoa, sideslip, &
                            gcr_max_projections, CFL

        use utils  , only : iflow_type, FLOW_INVISCID, FLOW_RANS, isolver_type, SOLVER_GCR, SOLVER_IMPLICIT, &
                            igrad_method, GRAD_LSQ, ilsq_stencil, LSQ_STENCIL_WVERTEX

        use solution_vars

        use turb

        implicit none

        ! initialize
        allocate( q(nq,ncells) )
        q = zero

        if (iflow_type > FLOW_INVISCID) allocate(mu(ncells))

        allocate( dtau(ncells) )
        allocate( wsn( ncells) )
        dtau = zero
        wsn = zero

        if ( accuracy_order > 1 .OR. iflow_type > FLOW_INVISCID) then
            allocate( ccgradq(ndim,nq,ncells) )
            if (igrad_method == GRAD_LSQ .and. ilsq_stencil == LSQ_STENCIL_WVERTEX) then
                allocate(  vgradq(ndim,nq,nnodes) )
            endif
        endif

        allocate( res(nq,ncells) )
        res = zero

        if ( isolver_type == SOLVER_IMPLICIT .or. isolver_type == SOLVER_GCR) then 
            allocate(solution_update(nq,ncells))
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

        CFL_used = CFL

        if (iflow_type >= FLOW_RANS) call allocate_rans

    end subroutine allocate_solution_vars

    subroutine compute_local_time_step_dtau

        use common                  , only : half
        use grid                    , only : ncells, cell
        use config                  , only : CFL, high_ar_correction
        use grid_statists           , only : cell_aspect_ratio
        use solution_vars           , only : wsn, dtau

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
    !  Input:  w = conservative variables (rho, u, v, w, p)
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

    use solution_vars , only : gamma
    
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
        use solution_vars , only : gamma, gmoinv, iu, iv, iw
        
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

    !********************************************************************************
    ! Compute rho from Q
    !
    ! ------------------------------------------------------------------------------
    ! This is just the first line of q2u 
    ! ------------------------------------------------------------------------------
    !
    ! Note: rho*E = p/(gamma-1) + rho*0.5*(u^2 + v^2 + w^2)
    !       rho   = p*gamma/T
    !********************************************************************************
    pure function q2rho(q_in) result(rho_out)

        use common, only : p2, half

        use solution_vars, only : gamma
        
        implicit none
    
        real(p2), dimension(5), intent(in) :: q_in ! input
        real(p2)                           :: rho_out !output
        
        ! rho    = p      *gamma / T
        rho_out = q_in(1)*gamma / q_in(5)
        
    end function q2rho

    pure function compute_primative_jacobian(qi) result(preconditioner)

        ! This function computes the Jacobian DU/DQ where U is the vector of conserved variables [rho rhoU rhoV rhoW rhoE] and W is
        ! the vector of primitive variables [p U V W T].  This is also written in a way that can allow implementation of Weiss-Smith
        ! Preconditioning

        use common , only : p2, half, zero, one

        use solution_vars , only : gmoinv, ip, iu, iv, iw, iT, gamma, gammamo

        implicit none

        real(p2), dimension(5),  intent(in) :: qi
        real(p2), dimension(5,5)            :: preconditioner
        
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
        preconditioner(5,:) = (/ theta*H-one ,  rho*qi(iu), rho*qi(iv), rho*qi(iw), rho_T*H + rho/(gamma-one)/)

    end function compute_primative_jacobian

end module solution

