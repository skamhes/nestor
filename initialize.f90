module initialize

    implicit none

    public :: set_initial_solution
    public :: init_jacobian
    public :: init_line_jacobian

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

        use config , only : line_implicit

        implicit none

        integer :: i, j, k
        integer :: c1, c2, length

        if (line_implicit) then
            call init_line_jacobian
            return
        endif

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

    ! The goal of the line Jacobian is to store it using the same structure (an array) as the point jacobian, that way it can be 
    ! handled with minimal additional control flow outside of this.
    !
    ! Each row of A will be split into two rows of non line elements and line elements.  The non line elements will be grouped first
    ! to update the off RHS of the array.  In additional to the normal R vector there is an additional row pointer that points to 
    ! the start of the off block and the line block.  The cells that don't belong to any line will be placed at the end of the 
    ! jacobian array and smoothed as normal.
    !
    !
    ! Storage pattern: A
    !
    !...              N-1                 N                     N+1               ...
    !Rline(2N)  =>|x                 |                    |  x        x       | R( 1) = 1
    !             |   x              |                    |     x             | R( 2) = 4
    !             |         x        |                    |   x               | R( 3) = 6
    !             |     x            |                    |    x              | R( 4) = 8
    !             |                  |                    |         x      x  | R( 5) = 10
    !             |         x  x     |                    |                   | R( 6) = 12
    !             |                  |                    |      x      x     | R( 7) = 14
    !             |          x  x    |                    |                   | R( 8) = 16
    !             |              x   |                    |       x           | R( 9) = 18
    !             |________________x_|____________________|x__________x_______| R(10) = 20
    !Rline(2N+1)=>|                  |d u                 |                   | R(11) = 23
    !             |                  |l d u               |                   | R(12) = 25
    !             |                  |  l d u             |                   | R(13) = 28
    !             |                  |    l d u           |                   | R(14) = 31
    !             |                  |      l d u         |                   | R(15) = 34
    !             |                  |        l d u       |                   | R(16) = 37
    !             |                  |          l d u     |                   | R(17) = 40
    !             |                  |            l d u   |                   | R(18) = 43
    !             |                  |              l d u |                   | R(19) = 46
    !             |                  |                l d |                   | R(20) = 49
    !Rline(2(N+1))=> ...                                                        R(21) = 51                                                                     
  
  
    subroutine init_line_jacobian

        use grid            , only : nfaces, face, cell, ncells

        use solution_vars        , only : nq, kth_nghbr_of_1, kth_nghbr_of_2, jac, diag_inv, C, R, nnz, kth_of_cell, Rline

        use sparse_common, only: insertion_sort_index

        use limplicit , only : nlines, lines

        implicit none

        integer :: i, j, k, jp
        integer :: cjm1, ck, cj
        integer :: c1, c2, length
        integer :: nrows

        integer, dimension(7)      :: nghbrs ! sorted scratch vector of cell neighbors and cell itself, 6 neighbors + 1
        integer, dimension(ncells) :: id_line

        id_line = 0

        ! Create kth_nghbr arrays
        if(.not.allocated(kth_nghbr_of_1)) allocate(kth_nghbr_of_1(nfaces))
        if(.not.allocated(kth_nghbr_of_2)) allocate(kth_nghbr_of_2(nfaces))
        if(.not.allocated(kth_of_cell)   ) allocate(kth_of_cell(ncells))

        allocate(Rline(2 * (nlines + 1)))
        nrows = 0

        do i = 1,nlines
            nrows = nrows + 2 * lines(i)%ncells
        end do

        allocate(R(ncells + 1 + nrows / 2))

        Rline(1) = 1
        R(1)     = 1
        nrows    = 1

        do i = 1,nlines
            if (i > 1) then
            end if

            ! line 2i + 2 of R 
            cjm1 = lines(i)%lcells(1)
            id_line(cjm1) = i
            R(nrows+1) = R(nrows) + cell(cjm1)%nnghbrs - 1 ! subtract the 1 neighbor going into the line block's 1st line
            nrows = nrows + 1
            
            do j = 3,lines(i)%ncells
                cjm1 = lines(i)%lcells(j-1)
                id_line(cjm1) = i
                R(nrows+1) = R(nrows) + cell(cjm1)%nnghbrs - 2 ! subtract the 2 neighbors going into the line block
                nrows = nrows + 1
            end do

            ! First line of the line block
            j = lines(i)%ncells + 1
            cjm1 = lines(i)%lcells(j-1)
            id_line(cjm1) = i
            R(nrows+1) = R(nrows) + cell(cjm1)%nnghbrs - 1 ! subtract the 1 neighbor going into the line block's last line
            nrows = nrows + 1
            
            Rline(2*i) = R(nrows)

            ! line 2 is also different
            R(nrows+1) = R(nrows) + 2 ! only two elemnts in the first line
            nrows = nrows + 1
            
            do j = 2,lines(i)%ncells-1
                R(nrows+1) = R(nrows) + 3 ! 3 elemnts in lines 2:(n-1)
                nrows = nrows + 1
            end do

            ! last line is also different
            R(nrows+1) = R(nrows) + 2 ! only two elemnts in the first line
            nrows = nrows + 1

            Rline(2 * i + 1) = R(nrows)
        end do

        ! Add the remaining cells to a standard jacobian block at the end
        do i = 1,ncells
            if (id_line(i) > 0) cycle
            R(nrows + 1) = R(nrows) + 1 + cell(i)%nnghbrs ! Start of row(nrows+1) = row(nrows) start point + 1 (diagonal term) + # of neighbors
            nrows = nrows + 1
        end do
        nnz = R(nrows) - 1 ! number of nonzero cells
        
        allocate(Jac(5,5,nnz))
        allocate(C(      nnz))
        
        allocate(diag_inv(5,5,ncells))
        
        do i = 1,nlines
            ! using 2 pointers to write the cell indices
            jp = Rline(2*i - 1)
            do k = 1,lines(i)%ncells
                ck = lines(i)%lcells(k)
                ! sort the index of the cell neighbors and i and stores them in C:
                length = cell(ck)%nnghbrs
                call insertion_sort_index( (/ cell(ck)%nghbr /) , nghbrs(1:length) ) 
                do j = 1,length
                    if(id_line( cell(ck)%nghbr( nghbrs(j) ) ) /= i) then ! not on the line
                        C(jp) = cell(ck)%nghbr(nghbrs(j))
                        jp = jp + 1
                    end if
                end do
            end do

            ! First row
            C(jp  ) = lines(i)%lcells(1)
            C(jp+1) = lines(i)%lcells(2)
            jp      = jp + 2
            ! Middle rows
            do k = 2,lines(i)%ncells - 1
                C(jp  ) = lines(i)%lcells(k-1)
                C(jp+1) = lines(i)%lcells(k  )
                C(jp+2) = lines(i)%lcells(k+1)
                jp      = jp + 3
            end do
            ! last row
            C(jp  ) = lines(i)%lcells(lines(i)%ncells - 1)
            C(jp+1) = lines(i)%lcells(lines(i)%ncells    )
            jp      = jp + 2
        end do        
        
        ! Create kth_nghbr arrays
        if(.not.allocated(kth_nghbr_of_1)) allocate(kth_nghbr_of_1(nfaces))
        if(.not.allocated(kth_nghbr_of_2)) allocate(kth_nghbr_of_2(nfaces))
        
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


    end subroutine init_line_jacobian
end module initialize