module forces

    implicit none

    public

    contains

    subroutine compute_forces

        use common      , only : p2, zero

        use config      , only : drag, lift, accuracy_order, turbulence_type, M_inf, Re_inf, sutherland_constant, reference_temp

        use utils       , only : ibc_type, BC_VISC_STRONG, iturb_type, TURB_INVISCID, ilsq_stencil, LSQ_STENCIL_WVERTEX

        use grid        , only : bound, nb, bc_type, cell

        use solution_vars, only : q, force_drag, force_lift, ndim, ccgradq, force_normalization, &
                                 vector_drag, vector_lift, vgradq

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
        real(p2), dimension(ndim,ndim) :: stress_tensor
        real(p2), dimension(ndim)      :: bface_shear
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
            if ( .NOT. (ibc_type(ib) == BC_VISC_STRONG .OR. trim(bc_type(ib)) == 'slip_wall' ) ) then
                ! Only one of these is integers for now.  We don't have a unique int for slip walls...
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

                if ( iflow_type > FLOW_INVISCID ) then
                    face_sides = bound(ib)%bfaces(1,i)

                    if (ilsq_stencil == LSQ_STENCIL_WVERTEX) then
                        bface_grad = zero
                        do k = 1,face_sides
                            nk = bound(ib)%bfaces(k + 1,i)
                            bface_grad = bface_grad + vgradq(:,2:4,nk)
                        end do
                        bface_grad = bface_grad / real(face_sides, p2)
                    else
                        bface_grad = ccgradq(1:3,2:4,ci)
                        ! This seems to be less accurate as the cell centered gradient is influenced by flow values farther from
                        ! the wall
                    endif
                    stress_tensor = compute_tau_wall(q(5,ci),bface_grad)
                    
                    bface_shear = matmul(stress_tensor,bface_normal)
                end if

                if ( lift ) then
                    pforce_lift = pforce_lift + face_pressure * bface_mag * dot_product(bface_normal,vector_lift)
                    if ( iflow_type > FLOW_INVISCID ) then 
                        vforce_lift = vforce_lift - bface_mag * dot_product(bface_shear, vector_lift)
                    end if
                endif
                if ( drag ) then
                    pforce_drag = pforce_drag + face_pressure * bface_mag * dot_product(bface_normal,vector_drag)
                    if ( iflow_type > FLOW_INVISCID ) then 
                        vforce_drag = vforce_drag - bface_mag * dot_product(bface_shear, vector_drag)
                    end if
                endif
            end do cell_loop
        end do boundary_loop

        if ( lift ) then
            force_lift = (vforce_lift + pforce_lift) * force_normalization
        endif
        if ( drag ) then
            force_drag = (vforce_drag + pforce_drag) * force_normalization
        endif


    end subroutine compute_forces

    subroutine output_forces

        use config , only : lift, drag

        use solution_vars , only : force_drag, force_lift

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

        use solution_vars , only : force_lift, force_drag

        implicit none

        write(*,*)
        write(*,*)
        if (lift) write(*,'(a,g0)') "  Lift force = ", force_lift
        if (drag) write(*,'(a,g0)') "  Drag force = ", force_drag
    end subroutine report_lift

    pure function compute_tau_wall(T,face_grad) result(tau)
        
        use common                  , only : p2, four_third, three_half, two_third, one

        use config                  , only : Pr, sutherland_constant, ideal_gas_constant, Re_inf, M_inf, reference_temp
        
        use solution_vars           , only : ndim, T_inf 
        implicit none

        real(p2),                 intent(in) :: T              ! Temperature at the attached cell
        real(p2), dimension(:,:), intent(in) :: face_grad      ! Face grad (mean of face node gradients)
        real(p2), dimension(ndim,ndim)       :: tau            ! Stress Tensor  
        
                ! Local Vars
        real(p2)                     :: mu
        real(p2)                     :: C0
        real(p2), dimension(ndim)    :: grad_u, grad_v, grad_w   !Interface gradients of velocities
        
        integer, parameter :: ix = 1
        integer, parameter :: iy = 2
        integer, parameter :: iz = 3

        C0= sutherland_constant/reference_temp
        mu =  M_inf/Re_inf * (one + C0/T_inf) / (T + C0/T_inf)*T**(three_half)

        grad_u = face_grad(:,ix)
        grad_v = face_grad(:,iy)
        grad_w = face_grad(:,iz)

        tau(ix,ix) =  mu*(four_third*grad_u(ix) - two_third*grad_v(iy) - two_third*grad_w(iz))
        tau(iy,iy) =  mu*(four_third*grad_v(iy) - two_third*grad_u(ix) - two_third*grad_w(iz))
        tau(iz,iz) =  mu*(four_third*grad_w(iz) - two_third*grad_u(ix) - two_third*grad_v(iy))
    
        tau(ix,iy) =  mu*(grad_u(iy) + grad_v(ix))
        tau(ix,iz) =  mu*(grad_u(iz) + grad_w(ix))
        tau(iy,iz) =  mu*(grad_v(iz) + grad_w(iy))
    
        tau(iy,ix) = tau(ix,iy)
        tau(iz,ix) = tau(ix,iz)
        tau(iz,iy) = tau(iy,iz)

    end function compute_tau_wall
    
end module forces