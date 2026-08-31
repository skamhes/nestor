module initialize

    implicit none

    public :: set_initial_solution
    public :: init_jacobian
    public :: init_line_jacobian

    contains

    subroutine set_initial_solution

        use common , only : p2

        use grid   , only : ncells

        use config , only : M_inf, aoa, sideslip, perturb_initial, random_perturb, lift, drag, area_reference, &
                            high_ar_correction, sutherland_constant, reference_temp, Re_inf, M_inf, restart, CFL, CFL_turb, &
                            accuracy_order

        use utils  , only : isolver_type, SOLVER_GCR, iflow_type, FLOW_INVISCID, FLOW_RANS

        use solution

        use gradient , only : init_gradients

        use grid_statists , only : init_ar_array, compute_aspect_ratio

        use viscosity , only : compute_viscosity

        use solution_vars , only : rho_inf, u_inf, v_inf, w_inf, p_inf, q, &
                                   CFL_used

        use turb , only : init_turb, turb_var, turb_res

        use inout , only : read_restart_file

        use wall_distance , only : compute_wall_distance

        implicit none

        integer                 :: i
#ifdef __INTEL_COMPILER
        real(p2) :: rval
#endif
        real(p2), dimension(5)  :: q_init

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
#ifdef __GFORTRAN__
                    q(2:4,i) = q(2:4,i) * rand(0)
#else
                    call random_number(rval)
                    q(2:4,i) = q(2:4,i) * rval
#endif
                endif
            end do cell_loop
        
        endif
        
        if (high_ar_correction) then
            call init_ar_array
            call compute_aspect_ratio
        endif

        if (iflow_type >= FLOW_RANS) then
            call init_turb
        else 
            nullify(turb_var, turb_res)
        endif

        if (isolver_type == SOLVER_GCR) then
            if (restart) then
                CFL = CFL_used
                CFL_turb = CFL_used
            else
                ! I'm not sure why but the GCR seems to behave better if the initial iterations are performed with a smaller CFL
                ! this not a very restrictive limitation as successful iterations will cause the CFL to grow quickly.
                CFL = min(CFL,1.0_p2)
                CFL_turb = min(CFL_turb,1.0_p2)
            end if
        end if

        if (accuracy_order == 2 .OR. iflow_type > FLOW_INVISCID ) then
            call init_gradients
        endif    

        if (iflow_type >= FLOW_RANS) call compute_wall_distance
    end subroutine set_initial_solution

    
    
    subroutine init_jacobian

        use grid            , only : nfaces, face, cell, ncells

        use solution_vars   , only : kth_nghbr_of_1, kth_nghbr_of_2, jac, diag_inv, c, R, nnz, kth_of_cell

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
    !Rline(2N-1)=>|x                 |                    |  x        x       | R( 1) = 1
    !             |   x              |                    |     x             | R( 2) = 4
    !             |         x        |                    |   x               | R( 3) = 6
    !             |     x            |                    |    x              | R( 4) = 8
    !             |                  |                    |         x      x  | R( 5) = 10
    !             |         x  x     |                    |                   | R( 6) = 12
    !             |                  |                    |      x      x     | R( 7) = 14
    !             |          x  x    |                    |                   | R( 8) = 16
    !             |              x   |                    |       x           | R( 9) = 18
    !             |________________x_|____________________|x__________x_______| R(10) = 20
    !Rline(2N)  =>|                  |d u                 |                   | R(11) = 23
    !             |                  |l d u               |                   | R(12) = 25
    !             |                  |  l d u             |                   | R(13) = 28
    !             |                  |    l d u           |                   | R(14) = 31
    !             |                  |      l d u         |                   | R(15) = 34
    !             |                  |        l d u       |                   | R(16) = 37
    !             |                  |          l d u     |                   | R(17) = 40
    !             |                  |            l d u   |                   | R(18) = 43
    !             |                  |              l d u |                   | R(19) = 46
    !             |                  |                l d |                   | R(20) = 49
    !Rline(2N+1)  => ...                                                        R(21) = 51                                                                     
  
  
    subroutine init_line_jacobian

        use grid            , only : nfaces, face, cell, ncells

        use solution_vars        , only : kth_nghbr_of_1, kth_nghbr_of_2, jac, diag_inv, C, R, nnz, kth_of_cell, Rline, iRow

        use sparse_common, only: insertion_sort_index

        use limplicit , only : nlines, lines

        implicit none

       
        integer, dimension(ncells) :: c2row ! for line cells, it is the offset from the base row,
        
        integer :: i, j, k, jp
        integer :: cjm1, ck, ci
        integer :: r1, r2, br
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

        allocate(R(ncells + 1 + nrows / 2)) ! Total number of rows counting the line cells twice
        allocate(iRow(nrows+1:ncells + nrows / 2)) ! arbitrary array indices in Fortran is neat

        Rline(1) = 1
        R(1)     = 1
        nrows    = 1

        do i = 1,nlines

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
            
            Rline(2*i) = nrows

            ! line 1 is also different
            R(nrows+1) = R(nrows) + 2 ! only two elemnts in the first line
            nrows = nrows + 1
            
            do j = 2,lines(i)%ncells-1
                R(nrows+1) = R(nrows) + 3 ! 3 elemnts in lines 2:(n-1)
                nrows = nrows + 1
            end do

            ! last line is also different
            R(nrows+1) = R(nrows) + 2 ! only two elemnts in the first line
            nrows = nrows + 1

            Rline(2 * i + 1) = nrows
        end do

        ! Add the remaining cells to a standard jacobian block at the end
        do i = 1,ncells
            if (id_line(i) > 0) cycle
            R(nrows + 1) = R(nrows) + 1 + cell(i)%nnghbrs ! Start of row(nrows+1) = row(nrows) start point + 1 (diagonal term) + # of neighbors
            id_line(i) = -nrows 
            iRow(nrows) = i
            nrows = nrows + 1
        end do
        Rline(2 * (nlines + 1)) = nrows
        nnz = R(nrows) - 1 ! number of nonzero cells
        
        allocate(Jac(5,5,nnz))
        allocate(C(      nnz))
        
        allocate(diag_inv(5,5,ncells))
        
        ! Write the C vector
        jp = 1
        do i = 1,nlines
            ! using 2 pointers to write the cell indices
            do k = 1,lines(i)%ncells
                ck = lines(i)%lcells(k)
                
                ! set the conversion c2row offset pointer
                c2row(ck) = k-1

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
            kth_of_cell(lines(i)%lcells(1)) = jp
            jp      = jp + 2
            ! Middle rows
            do k = 2,lines(i)%ncells - 1
                C(jp  ) = lines(i)%lcells(k-1)
                C(jp+1) = lines(i)%lcells(k  )
                C(jp+2) = lines(i)%lcells(k+1)
                kth_of_cell(lines(i)%lcells(k)) = jp + 1
                jp      = jp + 3
            end do
            ! last row
            C(jp  ) = lines(i)%lcells(lines(i)%ncells - 1)
            C(jp+1) = lines(i)%lcells(lines(i)%ncells    )
            kth_of_cell(lines(i)%lcells( lines(i)%ncells )) = jp + 1
            jp      = jp + 2
        end do    
        
        ! Add the rest
        do i = 1,ncells
            if (id_line(i) > 0) cycle
            ci = -id_line(i)
            ! sort the index of the cell neighbors and i and stores them in C:
            call insertion_sort_index( (/ cell(i)%nghbr, i /) , C(R(ci) : (R(ci+1)-1)) ) 
            length = R(ci+1)-R(ci)
            do j = R(ci),(R(ci+1)-1)
                if (length == C(j)) then
                    C(j) = i
                    kth_of_cell(i) = j
                else
                    C(j) = cell(i)%nghbr(C(j))
                end if
            end do
        end do
        
        ! Create kth_nghbr arrays
        if(.not.allocated(kth_nghbr_of_1)) allocate(kth_nghbr_of_1(nfaces))
        if(.not.allocated(kth_nghbr_of_2)) allocate(kth_nghbr_of_2(nfaces))

        ! Define kth neighbor arrays
        face_nghbr_loop : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            if (id_line(c1) == id_line(c2)) then ! in the same row
                br = Rline(2*id_line(c1)) ! baserow
                r1 = br + c2row(c1)
                r2 = br + c2row(c2)
            else ! not on the same line
                ! c1
                if (id_line(c1) > 0) then ! offline neighbor
                    br = Rline(2*id_line(c1)-1) ! baserow
                    r1 = br + c2row(c1)
                else
                    r1 = -id_line(c1) ! point neighbor
                endif
                !c2
                if (id_line(c2) > 0) then ! offline neighbor
                    br = Rline(2*id_line(c2)-1) 
                    r2 = br + c2row(c2)
                else
                    r2 = -id_line(c2) ! point neighbor
                endif
            endif

            ! loop over c1 neighbors to find c2
            k = findloc( C(R(r1) : R(r1+1)-1), c2, dim=1 )
            kth_nghbr_of_1(i) = R(r1) + k - 1

            k = findloc( C(R(r2) : R(r2+1)-1), c1, dim=1 )
            kth_nghbr_of_2(i) = R(r2) + k - 1
        end do face_nghbr_loop

    end subroutine init_line_jacobian
end module initialize