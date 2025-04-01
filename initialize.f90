module initialize

    implicit none

    public :: set_initial_solution
    public :: init_jacobian

    contains

    subroutine set_initial_solution

        use common , only : p2, one, pi, two

        use grid   , only : ncells

        use config , only : M_inf, aoa, sideslip, perturb_initial, random_perturb, lift, drag, area_reference, &
                            high_ar_correction, sutherland_constant, reference_temp, Re_inf, M_inf

        use utils  , only : isolver_type, SOLVER_GCR, SOLVER_IMPLICIT, iflow_type, FLOW_INVISCID, FLOW_RANS

        use solution

        use grid_statists , only : init_ar_array, compute_aspect_ratio

        use viscosity , only : C0, mu_norm, compute_viscosity

        use solution_vars , only : force_normalization, rho_inf, u_inf, v_inf, w_inf, p_inf, gamma, q, T_inf, mu_inf, mu

        use turb , only : init_turb

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
        
        if (isolver_type == SOLVER_IMPLICIT .OR. isolver_type == SOLVER_GCR ) call init_jacobian
        
        force_normalization = two / ( rho_inf * area_reference *  M_inf**2 )

        if (high_ar_correction) then
            call init_ar_array
            call compute_aspect_ratio
        endif
        
        if (iflow_type > FLOW_INVISCID) then
            C0 = (sutherland_constant/reference_temp) / T_inf
            mu_norm = M_inf/Re_inf
            mu_inf = compute_viscosity(T_inf)
            mu(:) = mu_inf
        end if

        if (iflow_type >= FLOW_RANS) call init_turb

    end subroutine set_initial_solution

    
    
    subroutine init_jacobian

        use grid            , only : nfaces, face, cell, ncells

        use solution_vars        , only : nq, jacobian_type, kth_nghbr_of_1, kth_nghbr_of_2, jac

        implicit none

        integer :: i, k
        integer :: c1, c2

        ! Create kth_nghbr arrays
        if(.not.allocated(kth_nghbr_of_1) )allocate(kth_nghbr_of_1(nfaces))
        if(.not.allocated(kth_nghbr_of_2) )allocate(kth_nghbr_of_2(nfaces))
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