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

        use grid            , only : ncells, nnodes, bc_type

        use least_squares   , only : lsq

        use solution        , only : ccgradq, vgradq, nq

        implicit none

        integer  :: i, ib, ivar
        integer  :: unknowns
        real(p2) :: wi

        vertex_loop : do i = 1, nnodes
            var_loop : do ivar = 1, nq
                if (lsq(i)%btype /= 0 ) then ! boundary vertex
                    ib = lsq(i)%btype
                    call boundary_value(bc_type(ib),ivar, unknowns, wi)
                endif
            end do var_loop    
        end do vertex_loop



        ! cell_loop : do i = 1,ncells
        !     var_loop : do ivar = 1,5
        !         wi = w(ivar,i)
        !         if (ivar == 5) wi = wi - p_inf + gauge_pressure
        !         nghbr_loop : do k = 1,cclsq(i)%nnghbrs_lsq
        !             nghbr_cell = cclsq(i)%nghbr_lsq(k)
        !             wk = w(ivar,nghbr_cell)
        !             if (ivar == 5) wk = wk - p_inf + gauge_pressure
        !             gradw(1,ivar,i) = gradw(1,ivar,i) + cclsq(i)%cx(k)*(wk-wi)
        !             gradw(2,ivar,i) = gradw(2,ivar,i) + cclsq(i)%cy(k)*(wk-wi)
        !             gradw(3,ivar,i) = gradw(3,ivar,i) + cclsq(i)%cz(k)*(wk-wi)
        !         end do nghbr_loop
        !     end do var_loop
        ! end do cell_loop

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