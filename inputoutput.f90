module inout

    implicit none

    public


    
    
    public :: write_tecplot_file_b
    ! public :: write_tecplot_file_v ! TODO

    contains 

    subroutine write_tecplot_file_b
        use common          , only : p2, zero
                                     
        use grid            , only : nnodes, x, y, z, &
                                     nb, &
                                     cell, &
                                     bound, bgrid_type

        use solution_vars   , only : q, nq, gamma

        use config          , only : project_name, io_path

        use lowlevel        , only : my_alloc_int_ptr

        use files
                                     
        implicit none 

        type bnode_type
                integer                               :: nbnodes 
                integer, dimension(:),  pointer       :: bnodes  
        end type bnode_type
        type(bnode_type), dimension(nb)               :: bnode_data
            
        type(bgrid_type), dimension(:), pointer       :: bound_export
        integer :: i, os
        
        integer                           :: j, k, ib, nk, j_count
        integer                           :: bcell_i, candidate_node
        real(p2), dimension(:,:), pointer :: qn
        real(p2), dimension(:), pointer   :: Mn, rhon
        integer , dimension(:  ), pointer :: nc
        logical                           :: already_added
        ! integer, dimension(:), pointer    :: nbnodes

        allocate(qn(nq, nnodes))
        allocate(nc(    nnodes))
        
        allocate(Mn(           nnodes))
        allocate(rhon(         nnodes))

        
        nc = 0
        qn = zero

        bound_loop : do ib = 1,nb
            bface_loop : do i = 1,bound(ib)%nbfaces
                bcell_i = bound(ib)%bcell(i)
                do k = 1,cell(bcell_i)%nvtx
                    qn(:,cell(bcell_i)%vtx(k)) = qn(:,cell(bcell_i)%vtx(k))  + q(:,bcell_i)
                    nc(cell(bcell_i)%vtx(k))   = nc(cell(bcell_i)%vtx(k)) + 1
                end do
            end do bface_loop
        end do bound_loop

        do j = 1,nnodes
            qn(:,j) = qn(:,j) / nc(j) ! copmute an average
            Mn(j) = sqrt(qn(2,j)**2 + qn(3,j)**2 + qn(4,j)**2)  ! mach number (wrt free stream a)
            ! rho    = p      * gamma / T
            rhon(j) = qn(1,j) * gamma / qn(5,j)
        end do

        allocate(bound_export(nb))
        boundary_loop : do ib = 1,nb
            bnode_data(ib)%nbnodes = zero
            allocate(bnode_data(ib)%bnodes(1))
            bnode_data(ib)%bnodes(1) = zero
            bound_export(ib)%nbfaces = bound(ib)%nbfaces
            allocate( bound_export(ib)%bfaces( 5,bound(ib)%nbfaces ) )  
            bound_export(ib)%bfaces = zero
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                bound_export(ib)%bfaces(1,i) = bound(ib)%bfaces(1,i)
                bface_vertex_loop : do k = 2,(bound(ib)%bfaces(1,i) + 1) ! loop through number of vertices for face
                    if (bnode_data(ib)%nbnodes == 0 ) then
                        bnode_data(ib)%nbnodes = bnode_data(ib)%nbnodes + 1
                        bnode_data(ib)%bnodes  = bound(ib)%bfaces(k,i)
                        bound_export(ib)%bfaces(k,i) = 1
                        cycle bface_vertex_loop 
                    end if 
                    candidate_node = bound(ib)%bfaces(k,i)
                    already_added = .false.
                    bnodes_loop : do nk = 1,bnode_data(ib)%nbnodes
                        if (candidate_node == bnode_data(ib)%bnodes(nk)) then
                            already_added = .true.
                            bound_export(ib)%bfaces(k,i) = nk
                            exit bnodes_loop  
                        end if 
                    end do bnodes_loop 
                    if (.not.already_added) then
                        bnode_data(ib)%nbnodes = bnode_data(ib)%nbnodes + 1
                        call my_alloc_int_ptr(bnode_data(ib)%bnodes,bnode_data(ib)%nbnodes)
                        bnode_data(ib)%bnodes(bnode_data(ib)%nbnodes) = candidate_node
                        bound_export(ib)%bfaces(k,i) = bnode_data(ib)%nbnodes
                    end if
                end do bface_vertex_loop 
            end do bfaces_loop
        end do boundary_loop

        write(*,*)
        write(*,*) "-------------------------------------------------------"
        write(*,*) ' Writing Tecplot file = ', trim(io_path)//trim(filename_tecplot_b)
        write(*,*)
    
        !Open the output file.
        open(unit=8, file=trim(io_path)//trim(filename_tecplot_b), status="unknown", iostat=os)   

        !---------------------------------------------------------------------------

        !(0)Header information

        write(8,*) 'TITLE = "GRID"'
        write(8,*) 'VARIABLES = "x","y","z","p","u","v","w","T","rho","M"'

        do ib = 1,nb
            write(8,*) 'ZONE T = "',ib,'"  n=', bnode_data(ib)%nbnodes, &
                            ',e=', bound(ib)%nbfaces,' , zonetype=fequadrilateral, datapacking=point'
            do j_count = 1,bnode_data(ib)%nbnodes
                j = bnode_data(ib)%bnodes(j_count)
                write(8,'(10es25.15)') x(j), y(j), z(j), qn(1,j), qn(2,j), qn(3,j), qn(4,j), qn(5,j), rhon(j), Mn(j)
            end do
            ! Loop through faces
            do i = 1,bound_export(ib)%nbfaces
                if (bound_export(ib)%bfaces(1,i) == 3) then ! write tri as a degenerate quad
                    write(8,'(4i10)') bound_export(ib)%bfaces(2,i), bound_export(ib)%bfaces(3,i), & 
                                        bound_export(ib)%bfaces(4,i), bound_export(ib)%bfaces(4,i)
                else if (bound_export(ib)%bfaces(1,i) == 4) then ! write quad
                    write(8,'(4i10)') bound_export(ib)%bfaces(2,i), bound_export(ib)%bfaces(3,i), & 
                    bound_export(ib)%bfaces(4,i), bound_export(ib)%bfaces(5,i)
                end if

            end do
        end do
        
        close(8)
        write(*,*)
        write(*,*) ' End of Writing Tecplot file = ', trim(filename_tecplot_b)
        write(*,*) "-------------------------------------------------------"
        write(*,*)
    end subroutine write_tecplot_file_b

    subroutine residual_status_header

        use config , only : lift, drag

        use utils  , only : isolver_type, SOLVER_GCR, SOLVER_IMPLICIT, iflow_type, iturb_model, FLOW_RANS, TURB_SA

        write(*,"(A)",advance="no") " Iteration   continuity   x-momemtum   y-momentum   z-momentum       energy"
        if (iflow_type == FLOW_RANS .and. iturb_model == TURB_SA) then
            write(*,"(A)",advance="no") "          nut"
        endif
        write(*,"(A)",advance="no") "      max-res"

        !Forces
        if (lift) write(*,"(A)",advance="no") "         lift"
        if (drag) write(*,"(A)",advance="no") "         drag"

        ! Split
        write(*,"(A)",advance="no") "   |"

        select case(isolver_type)
        case(SOLVER_IMPLICIT)
            write(*,*) " sweeps     reduction       time          CFL"
        case(SOLVER_GCR)
            write(*,*) " projs.     reduction       time          CFL"
        case default
            write(*,*) "   time     CFL"
        end select
    end subroutine residual_status_header

    subroutine print_residual_status(i_iteration, minutes, seconds)


        use config   , only : CFL, lift, drag

        use solution_vars , only : res_norm, res_norm_initial, lrelax_roc, lrelax_sweeps_actual, force_drag, force_lift, &
                              n_projections, nl_reduction, CFL_used

        use utils , only : isolver_type, SOLVER_IMPLICIT, SOLVER_GCR, iflow_type, iturb_model, FLOW_RANS, TURB_SA

        use turb  , only : turb_res_norm

        implicit none

        integer, intent(in) :: i_iteration
        integer, intent(in) :: minutes, seconds

        ! Formatting strings
        character(13) :: residual_format = '(i10,5es13.3)'
        character(8)  :: float_format    = '(es13.3)'
        character(45) :: implicit_format = '(i7,es14.1,i8.2,a,i2.2,es13.3)'
        character(35) :: explicit_format = '(i10.2,a,i2.2,es13.3)'


        ! Residuals
        write(*,residual_format,advance="no") i_iteration, res_norm(:)
        if (iflow_type == FLOW_RANS .and. iturb_model == TURB_SA) then
            write(*,float_format,advance="no") turb_res_norm
        endif
        write(*,float_format,advance="no") maxval(res_norm(:)/res_norm_initial(:))
        ! Forces
        if (lift) write(*,float_format,advance="no") force_lift
        if (drag) write(*,float_format,advance="no") force_drag

        write(*,"(A)",advance="no") "   | "



        ! Print out residual
        if ( isolver_type == SOLVER_IMPLICIT ) then
            write(*,implicit_format) lrelax_sweeps_actual, lrelax_roc, &
                                     minutes, ":", seconds, CFL
        elseif ( isolver_type == SOLVER_GCR ) then
            write(*,implicit_format) n_projections, nl_reduction, &
                                     minutes, ":", seconds, CFL_used
            ! note: the CFL is most likely wrong as it is updated for the next iteration before this line has had a chance to print.
            ! refactoring will take a little work so for now I'm leaving this note to remind me that it is wrong...
        else ! RK or Explicit
            write(*,explicit_format) minutes, ":", seconds, CFL
        endif

    end subroutine print_residual_status
    
end module inout