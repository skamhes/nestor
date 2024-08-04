module forces

    implicit none

    public

    contains

    subroutine compute_forces

        use common      , only : p2, zero

        use config      , only : drag, lift, accuracy_order

        use grid        , only : bound, nb, bc_type, cell

        use solution    , only : q, force_drag, force_lift, ndim, ccgradq, force_normalization, vector_drag, vector_lift

        implicit none

        integer                     :: ib, i, ci
        real(p2)                    :: cell_pressure
        real(p2), dimension(ndim)   :: pressure_grad
        real(p2)                    :: face_pressure
        real(p2)                    :: bface_mag
        real(p2), dimension(ndim)   :: bface_normal
        real(p2), dimension(ndim)   :: cell_center
        real(p2), dimension(ndim)   :: face_center
        real(p2), dimension(ndim)   :: cell_to_face_vect

        ! Initialize
        force_lift = zero
        force_drag = zero

        boundary_loop : do ib = 1,nb
            ! Only walls:
            if ( .NOT. (trim(bc_type(ib)) == 'no_slip_wall' .OR. trim(bc_type(ib)) == 'slip_wall' ) ) then
                cycle boundary_loop
            endif

            cell_loop : do i = 1,bound(ib)%nbfaces
                ci = bound(ib)%bcell(i)
                cell_pressure = q(1,ci)
                if ( accuracy_order == 2) then
                    pressure_grad(:) = ccgradq(:,1,ci)
                    cell_center = (/cell(ci)%xc,cell(ci)%yc,cell(ci)%zc/)
                    face_center = bound(ib)%bface_center(:,i)
                    cell_to_face_vect = face_center - cell_center
                    face_pressure = cell_pressure + dot_product(cell_to_face_vect,pressure_grad) ! extrapolate pressure to face
                else
                    face_pressure = cell_pressure
                endif
                bface_normal = bound(ib)%bface_nrml(:,i)
                bface_mag    = bound(ib)%bface_nrml_mag(i)
                if ( lift ) then
                    force_lift = force_lift + face_pressure * bface_mag * dot_product(bface_normal,vector_lift)
                endif
                if ( drag ) then
                    force_drag = force_drag + face_pressure * bface_mag * dot_product(bface_normal,vector_drag)
                endif
            end do cell_loop
        end do boundary_loop

        if ( lift ) force_lift = force_lift * force_normalization
        if ( drag ) force_drag = force_drag * force_normalization


    end subroutine compute_forces

    subroutine output_forces

        use config , only : lift, drag

        use solution , only : force_drag, force_lift

        ! character(11) :: format = '(a,es18.12)'
        character(11) :: format = '(g0)' !unlimited format specifier
        if (lift) then
            write(*,format) "Lift Force = ", force_lift
        endif
        if (drag) then
            write(*,format) "Drag Force = ", force_drag
        endif

    end subroutine output_forces
end module forces