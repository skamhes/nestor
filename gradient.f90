module gradient

    use common , only : p2

    implicit none 

    public init_gradients

    public compute_gradient_flow
    public compute_gradient_turb
    public set_ghost_values
    
    private


    contains

    subroutine init_gradients
        
        use least_squares , only : construct_lsq_stencil

        use config        , only : grad_method, lsq_stencil

        use solution_vars , only : ccgradq, vgradq

        use common        , only : zero

        use utils         , only : igrad_method, ilsq_stencil, GRAD_LSQ, LSQ_STENCIL_WVERTEX, iflow_type, FLOW_LAMINAR

        use turb            , only : ccgrad_turb_var, vgrad_turb_var


        implicit none
        
        select case(igrad_method)
        case(GRAD_LSQ)
            ! Build LSQ Stencil and coefficients
            call construct_lsq_stencil
            if (ilsq_stencil == LSQ_STENCIL_WVERTEX) then
                ! Gradient arrays are allocated in solution subroutine
                write(*,*)
                write(*,*) 'Initializing gradient vertex arrays.'
                vgradq = zero
                if (iflow_type > FLOW_LAMINAR) vgrad_turb_var = zero
            endif
        end select

        ! Initialize the cell centered array always
        write(*,*)
        write(*,*) 'Initializing gradient cell centered arrays.'
        write(*,*)
        
        ccgradq = zero
        if (iflow_type > FLOW_LAMINAR) ccgrad_turb_var = zero
        
    end subroutine init_gradients

    subroutine compute_gradient_flow(weight)
        ! This is a wrapper function to call the various gradient functions baced on the config.  Leaves the compute residual 
        ! subroutine a little cleaner...

        use common          , only : zero

        use config          , only : grad_method, lsq_stencil
        
        use utils           , only : igrad_method, ilsq_stencil, GRAD_LSQ, LSQ_STENCIL_WVERTEX, LSQ_STENCIL_NN

        use solution_vars   , only : ccgradq, vgradq

        implicit none

        integer, intent(in) :: weight ! this isn't used right now but eventually it should be

        integer :: dummy

        dummy = weight ! supress wunused At some point I may add weighted gradients (viscous terms?)
        
        ccgradq = zero

        select case(igrad_method)
        case(GRAD_LSQ)
            lsq : select case(ilsq_stencil)
            case(LSQ_STENCIL_WVERTEX) lsq
                vgradq = zero
                call compute_vgradient_flow
            case(LSQ_STENCIL_NN) lsq
                call compute_cgradient(weight)
            case default lsq
                write(*,*) 'Unsupported gradient methodstencil.'
                write(*,*) ' error in compute_gradients in gradient.f90. Stopping...'
                stop    
            end select lsq
        case default
            write(*,*) 'Unsupported gradient method.'
            write(*,*) ' error in compute_gradient_flow in gradient.f90. Stopping...'
            stop
        end select

        

    end subroutine compute_gradient_flow

    subroutine compute_vgradient_flow
        
        use common          , only : p2

        use grid            , only : ncells, nnodes, bound, cell

        use least_squares   , only : lsqv, INTERNAL

        use solution_vars   , only : ccgradq, vgradq, nq, q

        use bc_states       , only : get_right_state

        use utils           , only : ibc_type

        implicit none

        integer  :: i, ib, ivar, k, bound_int
        integer  :: attached_cell, attached_bface
        integer  :: unknowns
        real(p2) :: qi, qk          ! Node and attached cell values
        
        real(p2), dimension(5) :: qL, qcB              ! Values for attached bcell (for computing ghost cell values)
        real(p2), dimension(3) :: bface_nrml


        vertex_loop : do i = 1, nnodes
            var_loop : do ivar = 1, nq
                if (lsqv(i)%btype /= INTERNAL ) then ! boundary vertex
                    bound_int = lsqv(i)%btype
                    call boundary_value_flow(bound_int,ivar, unknowns, qi)
                else
                    unknowns = 4
                endif
                attach_loop : do k = 1,lsqv(i)%ncells_lsq
                    ! Get q at the neighboring cell (qk)
                    if (lsqv(i)%ib_lsq(k) == INTERNAL) then
                        attached_cell = lsqv(i)%cell_lsq(k)   
                        qk = q(ivar,attached_cell)
                    else ! BVERT
                        attached_bface = lsqv(i)%cell_lsq(k)
                        ib            = lsqv(i)%ib_lsq(k) ! this is the ib of the ghost cell (0 if internal)
                        attached_cell  = bound(ib)%bcell(attached_bface) 
                        qL = q(:,attached_cell)
                        bface_nrml = bound(ib)%bface_nrml(:,attached_bface)
                        ! This is somewhat redundant.  At some point I should improve it...
                        call get_right_state(qL,bface_nrml, ibc_type(ib),qcB)
                        qk = qcB(ivar)
                    endif
                    if ( unknowns == 3) then
                        vgradq(1,ivar,i) = vgradq(1,ivar,i) + lsqv(i)%cx3(k) * (qk - qi)
                        vgradq(2,ivar,i) = vgradq(2,ivar,i) + lsqv(i)%cy3(k) * (qk - qi)
                        vgradq(3,ivar,i) = vgradq(3,ivar,i) + lsqv(i)%cz3(k) * (qk - qi)
                    else
                        vgradq(1,ivar,i) = vgradq(1,ivar,i) + lsqv(i)%cx4(k) * qk
                        vgradq(2,ivar,i) = vgradq(2,ivar,i) + lsqv(i)%cy4(k) * qk
                        vgradq(3,ivar,i) = vgradq(3,ivar,i) + lsqv(i)%cz4(k) * qk
                    endif
                end do attach_loop

            end do var_loop    

            ! Add the vertex grad to each cell
            do k = 1,lsqv(i)%ncells_lsq
                if ( lsqv(i)%ib_lsq(k) == INTERNAL ) then
                    attached_cell = lsqv(i)%cell_lsq(k)
                    ccgradq(:,:,attached_cell) = ccgradq(:,:,attached_cell) + vgradq(:,:,i)
                endif
            enddo
        end do vertex_loop

        ! Devide cell gradient by number of attached vertices
        do i = 1,ncells
            ccgradq(:,:,i) = ccgradq(:,:,i) / real(cell(i)%nvtx, p2)
        end do

    end subroutine compute_vgradient_flow

    subroutine compute_cgradient(weight)

        use common , only : p2, ix, iy, iz

        use grid , only : nb, gcell, bound, ncells

        use solution_vars , only : q, ccgradq, nq, nlsq

        use utils , only : ibc_type

        use least_squares , only : lsqc
        
        implicit none

        integer, intent(in) :: weight

        integer :: ib, j, icell, kcell, jvar
        integer :: c1
        integer :: ck, ci

        real(p2), dimension(3) :: unit_face_normal
        real(p2), dimension(5) :: q1, qb
        real(p2), dimension(5) :: qk, qi, dq
        real(p2)               :: qk_j

        real(p2) :: dqx, dqy, dqz
        real(p2), dimension(3) :: dqf

        do icell=1,ncells
            qi = q(:,icell)
            do kcell = 1,lsqc(icell)%n_nnghbrs
                ck = lsqc(icell)%nghbr_lsq(kcell)
                qk = q(:,ck)
                dq(:) = qk - qi
                ! outer product
                do jvar = 1,5 ! Hard code for loop unrolling
                    ccgradq(:,jvar,icell) = ccgradq(:,jvar,icell) + lsqc(icell)%cf(:,kcell,weight) * dq(jvar)
                end do
                ! The loop unroller seems to be smarter than me doing this manually
                ! ccgradq(:,1,icell) = ccgradq(:,1,icell) + lsqc(icell)%cf(:,kcell,weight) * dq(1)
                ! ccgradq(:,2,icell) = ccgradq(:,2,icell) + lsqc(icell)%cf(:,kcell,weight) * dq(2)
                ! ccgradq(:,3,icell) = ccgradq(:,3,icell) + lsqc(icell)%cf(:,kcell,weight) * dq(3)
                ! ccgradq(:,4,icell) = ccgradq(:,4,icell) + lsqc(icell)%cf(:,kcell,weight) * dq(4)
                ! ccgradq(:,5,icell) = ccgradq(:,5,icell) + lsqc(icell)%cf(:,kcell,weight) * dq(5)
            end do
            do kcell = 1,lsqc(icell)%nbf
                ci = lsqc(icell)%gcells(1,kcell)
                ib = lsqc(icell)%gcells(2,kcell)
                qk = gcell(ib)%q(:,ci)
                dq(:) = qk - qi
                ! outer product
                do jvar = 1,5
                    ccgradq(:,jvar,icell) = ccgradq(:,jvar,icell) + lsqc(icell)%gcf(:,kcell,weight) * dq(jvar)
                end do
                ! ccgradq(:,1,icell) = ccgradq(:,1,icell) + lsqc(icell)%gcf(:,kcell,weight) * dq(1)
                ! ccgradq(:,2,icell) = ccgradq(:,2,icell) + lsqc(icell)%gcf(:,kcell,weight) * dq(2)
                ! ccgradq(:,3,icell) = ccgradq(:,3,icell) + lsqc(icell)%gcf(:,kcell,weight) * dq(3)
                ! ccgradq(:,4,icell) = ccgradq(:,4,icell) + lsqc(icell)%gcf(:,kcell,weight) * dq(4)
                ! ccgradq(:,5,icell) = ccgradq(:,5,icell) + lsqc(icell)%gcf(:,kcell,weight) * dq(5)
            end do
        end do


    end subroutine compute_cgradient

    subroutine boundary_value_flow(boundary_type, scalar, known, value)
        use common          , only : p2, zero

        use solution_vars   , only : p_inf, u_inf, v_inf, w_inf, T_inf

        use least_squares   , only : FREE_STREAM, SLIP_WALL, NO_SLIP_WALL, PRESSURE_OUTLET

        implicit none

        integer      ,              intent(in ) :: boundary_type
        integer      ,              intent(in ) :: scalar
        integer      ,              intent(out) :: known
        real(p2)     ,              intent(out) :: value
    
        
        select case(boundary_type) 
        case(FREE_STREAM)
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
        case(NO_SLIP_WALL)
            if (scalar >= 2 .AND. scalar <= 4) then
                value = zero
                known = 3
            else
                known = 4
            endif
        case(SLIP_WALL)
            ! For ghost values that depend on internal flow values, for now we will leave the vertex
            ! value as unknown
            known = 4
        case(PRESSURE_OUTLET)
            known = 4
        case default
            write(*,*) "Boundary condition #",boundary_type,"  not implemented."
            stop
    end select
        
    end subroutine boundary_value_flow

    subroutine compute_gradient_turb(weight)
        ! This is a wrapper function to call the various gradient functions baced on the config.  Leaves the compute residual 
        ! subroutine a little cleaner...

        use common          , only : zero

        use config          , only : grad_method, lsq_stencil
        
        use utils           , only : igrad_method, ilsq_stencil, GRAD_LSQ, LSQ_STENCIL_WVERTEX

        use turb            , only : ccgrad_turb_var, vgrad_turb_var

        implicit none

        integer, intent(in) :: weight ! this isn't used right now but eventually it should be

        integer :: dummy

        dummy = weight ! supress wunused At some point I may add weighted gradients (viscous terms?)
        
        ccgrad_turb_var = zero

        select case(igrad_method)
        case(GRAD_LSQ)
            if (ilsq_stencil == LSQ_STENCIL_WVERTEX) then
                vgrad_turb_var = zero

                call compute_vgradient_turb
            end if
        case default
            write(*,*) 'Unsupported gradient method.'
            write(*,*) ' error in compute_gradient_turb in gradient.f90. Stopping...'
            stop
        end select

    end subroutine compute_gradient_turb

    subroutine compute_vgradient_turb

        use grid            , only : ncells, nnodes, bound, cell

        use least_squares   , only : lsqv, INTERNAL

        use bc_states       , only : get_right_state

        use utils           , only : ibc_type

        use turb            , only : nturb, ccgrad_turb_var, vgrad_turb_var, turb_var

        use turb_bc         , only : sa_rhstate

        implicit none

        integer  :: i, ib, ivar, k, bound_int, icell
        integer  :: attached_cell, attached_bface
        integer  :: unknowns
        real(p2) :: ti, tk          ! Node and attached cell values
        
        real(p2) :: tL, tcB              ! Values for attached bcell (for computing ghost cell values)
        

        var_loop : do ivar = 1, nturb
            vertex_loop : do i = 1, nnodes
                if (lsqv(i)%btype /= INTERNAL ) then ! boundary vertex
                    bound_int = lsqv(i)%btype
                    call boundary_value_turb(bound_int, unknowns, ti)
                else
                    unknowns = 4
                endif
                attach_loop : do k = 1,lsqv(i)%ncells_lsq
                    ! Get q at the neighboring cell (tk)
                    if (lsqv(i)%ib_lsq(k) == INTERNAL) then
                        attached_cell = lsqv(i)%cell_lsq(k)   
                        tk = turb_var(attached_cell,ivar)
                    else ! BVERT
                        attached_bface = lsqv(i)%cell_lsq(k)
                        ib            = lsqv(i)%ib_lsq(k) ! this is the ib of the ghost cell (0 if internal)
                        attached_cell  = bound(ib)%bcell(attached_bface) 
                        tL = turb_var(attached_cell,ivar)
                        ! This is somewhat redundant.  At some point I should improve it...
                        call sa_rhstate(tL, ibc_type(ib), tcB)
                        tk = tcB
                    endif
                    if ( unknowns == 3) then
                        vgrad_turb_var(1,i,ivar) = vgrad_turb_var(1,i,ivar) + lsqv(i)%cx3(k) * (tk - ti)
                        vgrad_turb_var(2,i,ivar) = vgrad_turb_var(2,i,ivar) + lsqv(i)%cy3(k) * (tk - ti)
                        vgrad_turb_var(3,i,ivar) = vgrad_turb_var(3,i,ivar) + lsqv(i)%cz3(k) * (tk - ti)
                    else
                        vgrad_turb_var(1,i,ivar) = vgrad_turb_var(1,i,ivar) + lsqv(i)%cx4(k) * tk
                        vgrad_turb_var(2,i,ivar) = vgrad_turb_var(2,i,ivar) + lsqv(i)%cy4(k) * tk
                        vgrad_turb_var(3,i,ivar) = vgrad_turb_var(3,i,ivar) + lsqv(i)%cz4(k) * tk
                    endif
                end do attach_loop
            
                ! Add the vertex grad to each cell
                do k = 1,lsqv(i)%ncells_lsq
                    if ( lsqv(i)%ib_lsq(k) == INTERNAL ) then
                        attached_cell = lsqv(i)%cell_lsq(k)
                        ccgrad_turb_var(:,attached_cell, ivar) = ccgrad_turb_var(:,attached_cell,ivar) + vgrad_turb_var(:,i,ivar)
                    endif
                enddo
            end do vertex_loop

            
        end do var_loop    

        ! Devide cell gradient by number of attached vertices
        var_loop2 : do ivar = 1, nturb
            do icell = 1,ncells
                ccgrad_turb_var(:,icell,ivar) = ccgrad_turb_var(:,icell,ivar) / real(cell(icell)%nvtx, p2)
            end do
        end do var_loop2

    end subroutine compute_vgradient_turb

    subroutine boundary_value_turb(boundary_type, unknown, bval)
        use common          , only : p2, zero

        use least_squares   , only : FREE_STREAM, SLIP_WALL, NO_SLIP_WALL, PRESSURE_OUTLET

        use turb            , only : nut_inf

        use utils           , only : iturb_model, TURB_SA

        implicit none

        integer      ,              intent(in ) :: boundary_type
        integer      ,              intent(out) :: unknown
        real(p2)     ,              intent(out) :: bval
    
        if (iturb_model == TURB_SA) then
            select case(boundary_type) 
            case(FREE_STREAM)

            case(NO_SLIP_WALL)
                    bval = zero
            case(SLIP_WALL)
                ! For ghost values that depend on internal flow values, for now we will leave the vertex
                ! value as unknown
                unknown = 4
            case(PRESSURE_OUTLET)
                unknown = 3
                bval = nut_inf
            case default
                write(*,*) "Boundary condition #",boundary_type,"  not implemented."
                stop
            end select
        endif
    end subroutine boundary_value_turb

    subroutine set_ghost_values

        use solution_vars , only : q

        use bc_states , only : get_right_state
                
        use grid     , only : nb, bound, gcell

        use common   , only : p2

        use utils    , only : ibc_type

        implicit none

        integer  :: c1
        integer  :: ib, j
        
        real(p2), dimension(5) :: q1, qb
        real(p2), dimension(3) :: unit_face_normal

         ! First update the ghost cell values
        do ib = 1,nb
            do j=1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(j)
                unit_face_normal = bound(ib)%bface_nrml(:,j)
                q1 = q(:,c1)
                call get_right_state(q1, unit_face_normal, ibc_type(ib), qb)
                gcell(ib)%q(:,j) = qb
            end do
        end do

    end subroutine set_ghost_values

end module gradient