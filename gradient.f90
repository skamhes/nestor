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

        use utils         , only : igrad_method, ilsq_stencil, GRAD_LSQ, LSQ_STENCIL_WVERTEX

        !use grid          , only : nnodes, ncells


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
            endif
        end select

        ! Initialize the cell centered array always
        write(*,*)
        write(*,*) 'Initializing gradient cell centered arrays.'
        write(*,*)
        
        ccgradq = zero

    end subroutine init_gradients

    subroutine compute_gradient(weight)
        ! This is a wrapper function to call the various gradient functions baced on the config.  Leaves the compute residual 
        ! subroutine a little cleaner...

        use common          , only : p2, zero

        use config          , only : grad_method, lsq_stencil
        
        use utils           , only : igrad_method, ilsq_stencil, GRAD_LSQ, LSQ_STENCIL_WVERTEX, LSQ_STENCIL_NN

        use solution        , only : ccgradq, vgradq

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
                call compute_vgradient
            case(LSQ_STENCIL_NN) lsq
                call compute_cgradient
            case default lsq
                write(*,*) 'Unsupported gradient methodstencil.'
                write(*,*) ' error in compute_gradients in gradient.f90. Stopping...'
                stop    
            end select lsq
        case default
            write(*,*) 'Unsupported gradient method.'
            write(*,*) ' error in compute_gradients in gradient.f90. Stopping...'
            stop
        end select

        

    end subroutine compute_gradient

    subroutine compute_vgradient
        
        use common          , only : p2, zero

        use grid            , only : ncells, nnodes, bound, cell, x, y, z

        use least_squares   , only : lsqv, INTERNAL

        use solution        , only : ccgradq, vgradq, nq, q

        use bc_states       , only : get_right_state

        use utils           , only : ibc_type

        implicit none

        integer  :: i, ib, ivar, k, bound_int, cib
        integer  :: attached_cell, attached_bface
        integer  :: unknowns
        real(p2) :: qi, qk          ! Node and attached cell values
        real(p2) :: dx,dy,dz,cgx,cgy,cgz
        
        real(p2), dimension(5) :: qL, qcB              ! Values for attached bcell (for computing ghost cell values)
        real(p2), dimension(3) :: bface_nrml

        real(p2), dimension(3,5) :: dummy1, dummy2


        vertex_loop : do i = 1, nnodes
            var_loop : do ivar = 1, nq
                if (lsqv(i)%btype /= INTERNAL ) then ! boundary vertex
                    bound_int = lsqv(i)%btype
                    call boundary_value(bound_int,ivar, unknowns, qi, (/x(i),y(i),z(i)/))
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
                        cib = bound(ib)%bcell(attached_bface)
                        dx = bound(ib)%bface_center(1,attached_bface) - cell(cib)%xc
                        dy = bound(ib)%bface_center(2,attached_bface) - cell(cib)%yc
                        dz = bound(ib)%bface_center(3,attached_bface) - cell(cib)%zc
                        cgx = bound(ib)%bface_center(1,attached_bface) + dx
                        cgy = bound(ib)%bface_center(2,attached_bface) + dy
                        cgz = bound(ib)%bface_center(3,attached_bface) + dz
                        ! This is somewhat redundant.  At some point I should improve it...
                        call get_right_state(qL,(/cgx,cgy,cgz/),bface_nrml, ibc_type(ib),dummy1,qcB,dummy2)
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

    end subroutine compute_vgradient

    subroutine compute_cgradient

        use common , only : p2, ix, iy, iz

        use bc_states , only : get_right_state

        use grid , only : nb, gcell, bound, ncells, cell

        use solution , only : q, ccgradq, nq

        use utils , only : ibc_type

        use least_squares , only : lsqc
        
        implicit none

        integer :: ib, j, icell, kcell, jvar
        integer :: c1
        integer :: ck, ci

        real(p2), dimension(3) :: unit_face_normal
        real(p2), dimension(5) :: q1, qb
        real(p2), dimension(5) :: qk, qi
        real(p2)               :: qk_j
        
        real(p2), dimension(3,5) :: dummy1, dummy2

        real(p2)                :: xc, yc, zc
        real(p2)                :: fxc, fyc, fzc
        real(p2)                :: dxc2, dyc2, dzc2
        real(p2)                :: xc2, yc2, zc2
        

        ! First update the ghost cell values
        do ib = 1,nb
            do j=1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(j)
                unit_face_normal = bound(ib)%bface_nrml(:,j)
                q1 = q(:,c1)
                xc   = cell(c1)%xc
                yc   = cell(c1)%yc
                zc   = cell(c1)%zc
                fxc  = bound(ib)%bface_center(1,j)
                fyc  = bound(ib)%bface_center(2,j)
                fzc  = bound(ib)%bface_center(3,j)
                dxc2 = fxc - xc
                dyc2 = fyc - yc
                dzc2 = fzc - zc
                xc2  = fxc + dxc2
                yc2  = fyc + dyc2
                zc2  = fzc + dzc2
                call get_right_state(q1, (/xc2,yc2,zc2/), unit_face_normal, ibc_type(ib), dummy1, qb, dummy2)
                gcell(ib)%q(:,j) = qb
            end do
        end do


        do icell=1,ncells
            ! loop over the vertex neighboes
            qi = q(:,icell)
            do jvar = 1,nq
                do kcell = 1,lsqc(icell)%n_nnghbrs
                    ck = lsqc(icell)%nghbr_lsq(kcell)
                    qk_j = q(jvar,ck)
                    ccgradq(ix,jvar,icell) = ccgradq(ix,jvar,icell) + lsqc(icell)%cx(kcell) * (qk_j - qi(jvar))
                    ccgradq(iy,jvar,icell) = ccgradq(iy,jvar,icell) + lsqc(icell)%cy(kcell) * (qk_j - qi(jvar))
                    ccgradq(iz,jvar,icell) = ccgradq(iz,jvar,icell) + lsqc(icell)%cz(kcell) * (qk_j - qi(jvar))
                end do
                do kcell = 1,lsqc(icell)%nbf
                    ci = lsqc(icell)%gcells(1,kcell)
                    ib = lsqc(icell)%gcells(2,kcell)
                    qk_j = gcell(ib)%q(jvar,ci)
                    ccgradq(ix,jvar,icell) = ccgradq(ix,jvar,icell) + lsqc(icell)%gcx(kcell) * (qk_j - qi(jvar))
                    ccgradq(iy,jvar,icell) = ccgradq(iy,jvar,icell) + lsqc(icell)%gcy(kcell) * (qk_j - qi(jvar))
                    ccgradq(iz,jvar,icell) = ccgradq(iz,jvar,icell) + lsqc(icell)%gcz(kcell) * (qk_j - qi(jvar))
                end do
            end do

        end do


    end subroutine compute_cgradient

    subroutine boundary_value(boundary_type, scalar, known, bvalue,cc)
        use common          , only : p2, zero

        use solution        , only : p_inf, u_inf, v_inf, w_inf, T_inf

        use least_squares   , only : FREE_STREAM, SLIP_WALL, NO_SLIP_WALL, PRESSURE_OUTLET, MMS_DIRICHLET

        use mms             , only : fMMS

        implicit none

        integer      ,              intent(in ) :: boundary_type
        integer      ,              intent(in ) :: scalar
        integer      ,              intent(out) :: known
        real(p2)     ,              intent(out) :: bvalue

        real(p2),dimension(3),intent(in)  :: cc
        real(p2),dimension(5) ::qtmp, dummy
    
        
        select case(boundary_type) 
        case(FREE_STREAM)
            known = 3
            select case(scalar)
            case(1)
                bvalue = p_inf
            case(2)
                bvalue = u_inf
            case(3)
                bvalue = v_inf
            case(4)
                bvalue = w_inf
            case(5)
                bvalue = T_inf
            end select
        case(NO_SLIP_WALL)
            if (scalar >= 2 .AND. scalar <= 4) then
                bvalue = zero
                known = 3
            else
                known = 4
            endif
        case(SLIP_WALL)
            ! For ghost bvalues that depend on internal flow bvalues, for now we will leave the vertex
            ! bvalue as unknown
            known = 4
        case(PRESSURE_OUTLET)
            known = 4
        case(MMS_DIRICHLET)
            known = 3
            call fMMS(cc(1),cc(2),cc(3),qtmp,dummy)
            bvalue = qtmp(scalar)
        case default
            write(*,*) "Boundary condition #",boundary_type,"  not implemented."
            stop
    end select
        
    end subroutine boundary_value
end module gradient