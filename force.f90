module forces

    implicit none

    public

    contains

    subroutine compute_forces

        use common      , only : p2, zero, one, three_half, two, half

        use config      , only : drag, lift, accuracy_order, turbulence_type, M_inf, Re_inf, sutherland_constant, reference_temp

        use grid        , only : bound, nb, bc_type, cell

        use solution    , only : q, force_drag, force_lift, ndim, ccgradq, force_normalization, &
                                 vector_drag, vector_lift, vgradq, T_inf

        implicit none

        real(p2)                       :: cell_pressure
        real(p2), dimension(ndim)      :: pressure_grad
        real(p2)                       :: face_pressure
        real(p2)                       :: bface_mag
        real(p2), dimension(ndim)      :: bface_normal
        real(p2), dimension(ndim)      :: cell_center
        real(p2), dimension(ndim)      :: face_center
        real(p2), dimension(ndim)      :: cell_to_face_vect
        real(p2), dimension(ndim,ndim) :: bface_grad
        real(p2), dimension(ndim)      :: b_grad_vel_mag
        real(p2), dimension(ndim,ndim) :: stress_tensor
        real(p2), dimension(ndim)      :: bface_shear
        real(p2)                       :: T, C0, mu
        real(p2)                       :: pforce_lift, vforce_lift
        real(p2)                       :: pforce_drag, vforce_drag
        integer                        :: face_sides

        integer                        :: ib, i, ci, k, nk

        ! Initialize
        pforce_lift = zero
        pforce_drag = zero
        vforce_lift = zero
        vforce_drag = zero

        bface_shear = zero

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

                if ( trim(turbulence_type) /= 'inviscid' ) then
                    face_sides = bound(ib)%bfaces(1,i)

                    bface_grad = zero
                    do k = 1,face_sides
                        nk = bound(ib)%bfaces(k + 1,i)
                        bface_grad = bface_grad + vgradq(:,2:4,nk)
                    end do
                    bface_grad = bface_grad / real(face_sides, p2)
                    
                    stress_tensor = compute_tau_wall(q(5,ci),bface_grad)
                    
                    bface_shear = matmul(stress_tensor,bface_normal)
                end if

                if ( lift ) then
                    pforce_lift = pforce_lift + face_pressure * bface_mag * dot_product(bface_normal,vector_lift)
                    if ( trim(turbulence_type) /= 'inviscid' ) then 
                        vforce_lift = vforce_lift - bface_mag * dot_product(bface_shear, vector_lift)
                    end if
                endif
                if ( drag ) then
                    pforce_drag = pforce_drag + face_pressure * bface_mag * dot_product(bface_normal,vector_drag)
                    if ( trim(turbulence_type) /= 'inviscid' ) then 
                        vforce_drag = vforce_drag - bface_mag * dot_product(bface_shear, vector_drag)
                    end if
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

    subroutine report_lift

        use config , only : lift, drag

        use solution , only : force_lift, force_drag

        implicit none

        write(*,*)
        write(*,*)
        if (lift) write(*,'(a,g0)') "  Lift force = ", force_lift
        if (drag) write(*,'(a,g0)') "  Drag force = ", force_drag
    end subroutine report_lift

    function compute_tau_wall(T,face_grad) result(tau)
        
        use common                  , only : p2, half, one, four_third, three_half, two_third

        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp
        
        use solution                , only : gammamo, nq, ndim, T_inf 
        implicit none

        real(p2),                     intent(in) :: T              ! Temperature at the attached cell
        real(p2), dimension(ndim,nq), intent(in) :: face_grad      ! Face grad (mean of face node gradients)
        real(p2), dimension(ndim,ndim)           :: tau            ! Stress Tensor  
        
                ! Local Vars
        real(p2)                     :: mu
        real(p2)                     :: C0
        real(p2), dimension(3)       :: grad_u, grad_v, grad_w   !Interface gradients of velocities
        real(p2), dimension(3)       :: grad_T
        
        integer, parameter :: ix = 1
        integer, parameter :: iy = 2
        integer, parameter :: iz = 3

        C0= sutherland_constant/reference_temp
        mu =  M_inf/Re_inf * (one + C0/T_inf) / (T + C0/T_inf)*T**(three_half)

        grad_u = face_grad(:,2)
        grad_v = face_grad(:,3)
        grad_w = face_grad(:,4)

        tau(ix,ix) =  mu*(four_third*grad_u(1) - two_third*grad_v(2) - two_third*grad_w(3))
        tau(iy,iy) =  mu*(four_third*grad_v(2) - two_third*grad_u(1) - two_third*grad_w(3))
        tau(iz,iz) =  mu*(four_third*grad_w(3) - two_third*grad_u(1) - two_third*grad_v(2))
    
        tau(ix,iy) =  mu*(grad_u(2) + grad_v(1))
        tau(ix,iz) =  mu*(grad_u(3) + grad_w(1))
        tau(iy,iz) =  mu*(grad_v(3) + grad_w(2))
    
        tau(iy,ix) = tau(ix,iy)
        tau(iz,ix) = tau(ix,iz)
        tau(iz,iy) = tau(iy,iz)

    end function compute_tau_wall
    
end module forces