module gradient

    use common , only : p2

    implicit none 

    public init_gradients

    public compute_gradient

    contains

    subroutine init_gradients
        
        use least_squares , only : construct_lsq_stencil

        use config        , only : grad_method, lsq_stencil

        use solution      , only : ccgradq, vgradq

        use common        , only : p2, zero

        !use grid          , only : nnodes, ncells


        implicit none
        
        if (trim(grad_method) == 'lsq') then
            ! Build LSQ Stencil and coefficients
            call construct_lsq_stencil
            if (trim(lsq_stencil) == 'w_vertex') then
                ! Gradient arrays are allocated in solution subroutine
                write(*,*)
                write(*,*) 'Initializing gradient vertex arrays.'
                vgradq = zero
            endif
        endif

        ! Initialize the cell centered array always
        write(*,*)
        write(*,*) 'Initializing gradient cell centered arrays.'
        ccgradq = zero

    end subroutine init_gradients

    subroutine compute_gradient

        use common          , only : p2, zero

        use config          , only : grad_method, lsq_stencil

        use grid            , only : ncells, nnodes

        use least_squares   , only : lsq

        use solution        , only : ccgradq, vgradq

        ccgradq = zero

        if (trim(grad_method) == 'lsq') then
            if (trim(lsq_stencil) == 'w_vertex') then
                vgradq = zero

                call compute_vgradient
            end if
        else
            write(*,*) 'Unsupported gradient method.'
        endif

        

    end subroutine compute_gradient

    subroutine compute_vgradient
        
        use common          , only : p2, zero

        use grid            , only : ncells, nnodes, bc_type, bound

        use least_squares   , only : lsq, INTERNAL

        use solution        , only : ccgradq, vgradq, nq, q

        use bc_states       , only : get_right_state

        implicit none

        integer  :: i, ib, ivar, k
        integer  :: attached_cell, attached_bface
        integer  :: unknowns
        real(p2) :: qi, qk          ! Node and attached cell values
        
        real(p2), dimension(5) :: qL, qcB              ! Values for attached bcell (for computing ghost cell values)
        real(p2), dimension(3) :: bface_nrml


        vertex_loop : do i = 1, nnodes
            var_loop : do ivar = 1, nq
                if (lsq(i)%btype /= INTERNAL ) then ! boundary vertex
                    ib = lsq(i)%btype
                    call boundary_value(bc_type(ib),ivar, unknowns, qi)
                endif
                attach_loop : do k = 1,lsq(i)%ncells_lsq
                    if (lsq(i)%btype == INTERNAL) then
                        attached_cell = lsq(i)%cell_lsq(k)   
                        qk = q(ivar,attached_cell)
                        ! Add value to gradien
                        vgradq(1,ivar,i) = vgradq(1,ivar,i) + lsq(i)%cx4(k) * qk ! qi is on the LHS 
                        vgradq(2,ivar,i) = vgradq(2,ivar,i) + lsq(i)%cy4(k) * qk
                        vgradq(3,ivar,i) = vgradq(3,ivar,i) + lsq(i)%cz4(k) * qk
                    else ! BVERT
                        attached_bface = lsq(i)%cell_lsq(k)
                        ib            = lsq(i)%ib_lsq(k) ! this is the face ib not the node ib above
                        attached_cell  = bound(ib)%bcell(attached_bface) 
                        qL = q(:,attached_cell)
                        bface_nrml = bound(ib)%bface_nrml(:,attached_bface)
                        ! This is somewhat redundant.  At some point I should improve it...
                        call get_right_state(qL,bface_nrml, bc_type(ib), qcB)
                        qk = qcB(ivar)
                        if ( unknowns == 3) then
                            vgradq(1,ivar,i) = vgradq(1,ivar,i) + lsq(i)%cx3(k) * (qk - qi)
                            vgradq(2,ivar,i) = vgradq(2,ivar,i) + lsq(i)%cy3(k) * (qk - qi)
                            vgradq(3,ivar,i) = vgradq(3,ivar,i) + lsq(i)%cz3(k) * (qk - qi)
                        else
                            vgradq(1,ivar,i) = vgradq(1,ivar,i) + lsq(i)%cx4(k) * qk
                            vgradq(2,ivar,i) = vgradq(2,ivar,i) + lsq(i)%cy4(k) * qk
                            vgradq(3,ivar,i) = vgradq(3,ivar,i) + lsq(i)%cz4(k) * qk
                        endif
                    endif
                end do attach_loop

            end do var_loop    

            ! Add the vertex grad to each cell
            do k = 1,lsq(i)%ncells_lsq
                if ( lsq(i)%ib_lsq(k) == INTERNAL ) then
                    attached_cell = lsq(i)%cell_lsq(k)
                    ccgradq(:,:,attached_cell) = ccgradq(:,:,attached_cell) + vgradq(:,:,i)
                endif
            enddo
        end do vertex_loop

        ! Devide cell gradient by number of attached vertices
        do i = 1,ncells
            ccgradq(:,:,i) = ccgradq(:,:,i) / real(cell(i)%nvtx, p2)
        end do

    end subroutine compute_vgradient

    subroutine boundary_value(boundary_type, scalar, known, value)
        use common      , only : p2, zero

        use solution    , only : p_inf, u_inf, v_inf, w_inf, T_inf

        use grid        , only : bc_type

        implicit none

        character(80),              intent(in ) :: boundary_type
        integer      ,              intent(in ) :: scalar
        integer      ,              intent(out) :: known
        real(p2)     ,              intent(out) :: value
    
        
        select case(trim(boundary_type)) 
        case('freestream')
            known = 3
            select case(scalar)
            case(1)
                value = p_inf
            case(2)
                value = u_inf
            case(3)
                value = v_inf
            case(4)
                value = w_inf
            case(5)
                value = T_inf
            end select
        case('slip_wall')
            if (scalar >= 2 .AND. scalar <= 4) then
                value = zero
                known = 3
            else
                known = 4
            endif
        case('no_slip_wall')
            if (scalar == 1) then
                value = p_inf
                known = 3
            else
                ! For ghost values that depend on internal flow values, for now we will leave the vertex
                ! value as unknown
                known = 4
            endif
        case('outflow_subsonic')
            known = 4
        case default
            write(*,*) "Boundary condition=",trim(boundary_type),"  not implemented."
            stop
    end select
        
    end subroutine boundary_value
end module gradient