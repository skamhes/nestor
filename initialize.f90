module initialize

    implicit none

    public :: set_initial_solution
    public :: init_jacobian

    contains

    subroutine set_initial_solution

        use common , only : p2, one, pi, two

        use grid   , only : ncells

        use config , only : M_inf, aoa, sideslip, perturb_initial, random_perturb, lift, drag, area_reference, &
                            high_ar_correction, sutherland_constant, reference_temp, Re_inf, M_inf, restart

        use utils  , only : isolver_type, SOLVER_GCR, SOLVER_IMPLICIT, iflow_type, FLOW_INVISCID, FLOW_RANS

        use solution

        use grid_statists , only : init_ar_array, compute_aspect_ratio

        use viscosity , only : compute_viscosity

        use solution_vars , only : force_normalization, rho_inf, u_inf, v_inf, w_inf, p_inf, gamma, q, T_inf, mu_inf, mre, C0

        use turb , only : init_turb, nturb, turb_var, turb_res

        use inout , only : read_restart_file

        implicit none

        integer                 :: i
        real(p2), dimension(5)  :: q_init

        ! Set the free stream values
        rho_inf = one
        u_inf = M_inf*cos(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        v_inf = M_inf*sin(sideslip*pi/180)
        w_inf = M_inf*sin(aoa*pi/180_p2)*cos(sideslip*pi/180_p2)
        p_inf = one/gamma

        if (restart) then ! annoyingly I have a bunch of allocations inside the initialization subroutine. I'll have to seperate them out...
            call read_restart_file
        else
            
            
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
        
        endif
        
        
        force_normalization = two / ( rho_inf * area_reference *  M_inf**2 )

        if (high_ar_correction) then
            call init_ar_array
            call compute_aspect_ratio
        endif
        
        ! if (iflow_type > FLOW_INVISCID) then
        !     C0 = (sutherland_constant/reference_temp) / T_inf
        !     mu_norm = M_inf/Re_inf
        !     mu_inf = compute_viscosity(T_inf)
        !     mu(:) = mu_inf
        ! end if

        C0 = sutherland_constant/reference_temp
        mre = M_inf / Re_inf
        mu_inf = compute_viscosity(T_inf)

        if (iflow_type >= FLOW_RANS) then
            call init_turb
        else 
            nullify(turb_var, turb_res)
            nturb = 0 ! gonna use this in the GCR to be a little clever
        endif

    end subroutine set_initial_solution

    
    
    subroutine init_jacobian

        use grid            , only : nfaces, face, cell, ncells

        use solution_vars        , only : nq, kth_nghbr_of_1, kth_nghbr_of_2, jac, diag_inv, c, R, nnz, kth_of_cell

        use sparse_common, only: insertion_sort_index


        implicit none

        integer :: i, j, k
        integer :: c1, c2, length

        ! Create kth_nghbr arrays
        if(.not.allocated(kth_nghbr_of_1)) allocate(kth_nghbr_of_1(nfaces))
        if(.not.allocated(kth_nghbr_of_2)) allocate(kth_nghbr_of_2(nfaces))
        if(.not.allocated(kth_of_cell)   ) allocate(kth_of_cell(ncells))

        ! Count the number of nonzero values
        allocate(R(ncells + 1))
        R(1) = 1
        do i = 2,ncells + 1
            R(i) = R(i-1) + 1 + cell(i-1)%nnghbrs ! Start of row(i) = row(i-1) start point + 1 (diagonal term) + # of neighbors
        end do
        nnz = R(ncells+1) - 1 ! number of nonzero cells

        allocate(Jac(5,5,nnz))
        allocate(C(      nnz))

        allocate(diag_inv(5,5,ncells))

        do i = 1,ncells
            ! sort the index of the cell neighbors and i and stores them in C:
            call insertion_sort_index( (/ cell(i)%nghbr, i /) , C(R(i) : (R(i+1)-1)) ) 
            length = R(i+1)-R(i)
            do j = R(i),(R(i+1)-1)
                if (length == C(j)) then
                    C(j) = i
                    kth_of_cell(i) = j
                else
                    C(j) = cell(i)%nghbr(C(j))
                end if
            end do
        end do        
        
        ! Define kth neighbor arrays
        face_nghbr_loop : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            ! loop over c1 neighbors to find c2
            k = findloc( C(R(c1) : R(c1+1)-1), c2, dim=1 )
            kth_nghbr_of_1(i) = R(c1) + k - 1

            k = findloc( C(R(c2) : R(c2+1)-1), c1, dim=1 )
            kth_nghbr_of_2(i) = R(c2) + k - 1

        end do face_nghbr_loop

    end subroutine init_jacobian
end module initialize