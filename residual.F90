module residual

    implicit none 
    
    public :: compute_residual

    contains 

    subroutine compute_residual

        use utils , only : iflow_type, FLOW_RANS

        use res_turb

        call compute_residual_flow

        if (iflow_type == FLOW_RANS) call compute_residual_turb

    end subroutine compute_residual

    subroutine compute_residual_flow

        use common          , only : p2, zero, one, half

        use config          , only : method_inv_flux, accuracy_order, use_limiter

        use utils           , only : iflow_type, FLOW_INVISCID, FLOW_RANS

        use grid            , only : ncells, cell,  &
                                     nfaces, face,  &
                                     nb,     bound, &
                                     face_nrml,     &
                                     face_nrml_mag, &
                                     face_centroid

        use utils           , only : ibc_type
        
        use solution_vars   , only : res, q, ccgradq, vgradq, wsn, phi, mu, iT

        use solution        , only : q2u

        use interface       , only : interface_flux

        use limiter         , only : compute_limiter_flow

        use bc_states       , only : get_right_state

        use turb_bc         , only : turb_rhstate

        use gradient        , only : compute_gradient_flow

        use viscosity       , only : compute_viscosity

        use viscous_flux    , only : visc_flux_boundary, visc_flux_internal

        use turb            , only : turb_var, calcmut, nturb

        implicit none

        ! Grid Vars
        integer                     :: c1, c2
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid

        ! Flow variables
        real(p2), dimension(5)      :: q1, q2
        real(p2)                    :: mu1, mu2, muf, mutf
        real(p2), dimension(nturb)  :: trbv1, trbv2
        real(p2), dimension(3,5)    :: gradq1, gradq2, gradqb
        real(p2), dimension(5)      :: num_flux
        real(p2), dimension(5)      :: qb
        real(p2)                    :: wave_speed
        real(p2)                    :: phi1, phi2

        ! Misc int/counters
        integer                     :: i
        integer                     :: j, ib
        integer                     :: face_sides
        integer                     :: k, nk

        
        !--------------------------------------------------------------------------------
        ! Initialize the residuals and wsn = the sum of (max_wave_speed)*(face length)).
        ! Note: wsn is required to define a time step.
        cell_loop1 :  do i = 1, ncells

            res(:,i) = zero
            wsn(i)   = zero
     
        end do cell_loop1

        !--------------------------------------------------------------------------------
        ! Compute gradients at cells.
        !
        if (iflow_type > FLOW_INVISCID) then
            do i = 1,ncells
                mu(i) = compute_viscosity(q(iT,i))
            end do
            call compute_gradient_flow(0) ! For now we are just using unweighted gradients
        elseif ( accuracy_order == 2 ) then
            call compute_gradient_flow(0)
        endif

        if (use_limiter) call compute_limiter_flow

        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        !--------------------------------------------------------------------------------
        ! Residual computation: interior faces

        !--------------------------------------------------------------------------------
        ! Flux computation across internal faces (to be accumulated in res(:))
        !
        !          v2=Left(2)
        !        o---o---------o       face(j,:) = [i,k,v2,v1]
        !       .    .          .
        !      .     .           .
        !     .      .normal      .
        !    .  Left .--->  Right  .
        !   .   c1   .       c2     .
        !  .         .               .
        ! o----------o----------------o
        !          v1=Right(1)
        !
        !
        ! 1. Extrapolate the solutions to the face-midpoint from centroids 1 and 2.
        ! 2. Compute the numerical flux.
        ! 3. Add it to the residual for 1, and subtract it from the residual for 2.
        !
        !--------------------------------------------------------------------------------
        mutf = 0 ! only need to set once if not used

        loop_faces : do i = 1,nfaces
            ! Left and right cell values
            c1 = face(1,i)
            c2 = face(2,i)
            q1 = q(1:5, c1)
            q2 = q(1:5, c2)
            if (accuracy_order == 2 ) then
                gradq1 = ccgradq(1:3,1:5,c1)
                gradq2 = ccgradq(1:3,1:5,c2)
            else
                gradq1 = zero
                gradq2 = zero
            endif
            ! Face normal
            unit_face_normal = face_nrml(1:3,i)
            ! Limiters
            if (use_limiter) then
                phi1 = phi(c1)
                phi2 = phi(c2)
            else 
                phi1 = one
                phi2 = one
            end if
            call interface_flux(          q1,       q2   , & !<- Left/right states
                                      gradq1,      gradq2, & !<- Left/right gradients
                                         unit_face_normal, & !<- unit face normal
                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                    cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                       face_centroid(1,i), &
                                       face_centroid(2,i), &
                                       face_centroid(3,i), & !<- face midpoint
                                        phi1,        phi2, & !<- Limiter functions
                                     num_flux, wave_speed  ) !<- Output
 
            res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)
            wsn(c1)   = wsn(c1) + wave_speed*face_nrml_mag(i)
            
            res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)
            wsn(c2)   = wsn(c2) + wave_speed*face_nrml_mag(i)

            if ( iflow_type == FLOW_INVISCID) cycle loop_faces

            mu1 = mu(c1)
            mu2 = mu(c2)
            muf = half * (mu1 + mu2) ! we do this here because we need it more than once
            if (iflow_type == FLOW_RANS) then
                trbv1 = turb_var(c1,:)
                trbv2 = turb_var(c2,:)  ! don't like the non-contiguous memory reads but we're gonna have to do it somewhere.  
                ! Maybe I'll change it in the future...
                mutf = calcmut(q1,q2,muf,trbv1,trbv2)
                ! no elseif needed, we set this to zero before the loop.
            end if

            ! Viscous flux
            call visc_flux_internal(q1,q2,muf,mutf,gradq1,gradq2,unit_face_normal,  &
                                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                                    cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, &
                                                                   num_flux)

            res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)

            res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)

#ifdef NANCHECK
            if (any(isnan(res(:,c1)))) then 
                write (*,*) "nan value present - press [Enter] to continue"
                read(unit=*,fmt=*)
            end if
#endif

        end do loop_faces

        boundary_loop : do ib = 1,nb
            bface_loop : do j = 1,bound(ib)%nbfaces
                bface_centroid = bound(ib)%bface_center(:,j)
                if (use_limiter) then
                    phi1 = phi(c1)
                    phi2 = one
                else 
                    phi1 = one
                    phi2 = one
                end if

                c1 = bound(ib)%bcell(j)
                                
                unit_face_normal = bound(ib)%bface_nrml(:,j)
                
                gradq2 = zero ! won't matter since boundary cell center will be at face center

                q1 = q(:,c1)
                
                ! Get the right hand state (weak BC!)
                call get_right_state(q1, unit_face_normal, ibc_type(ib), qb)
                if ( accuracy_order == 2 ) then
                    gradq1 = ccgradq(1:3,1:5,c1)
                else
                    gradq1 = zero
                endif
                
                call interface_flux(          q1,      qb, & !<- Left/right states
                                         gradq1,   gradq2, & !<- Left/right gradients
                                         unit_face_normal, & !<- unit face normal
                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                                        bface_centroid(1), &
                                        bface_centroid(2), &
                                        bface_centroid(3), & !<- Right cell centroid
                                        bface_centroid(1), &
                                        bface_centroid(2), &
                                        bface_centroid(3), & !<- boundary ghost cell "center"
                                        phi1,        phi2, & !<- Limiter functions
                
                                        num_flux, wave_speed  )

                res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)
                wsn(c1)   = wsn(c1) + wave_speed * bound(ib)%bface_nrml_mag(j)

#ifdef NANCHECK
                if (any(isnan(res(:,c1)))) then 
                    write (*,*) "nan value present - press [Enter] to continue"
                    read(unit=*,fmt=*)
                end if
#endif
                if ( iflow_type == FLOW_INVISCID ) cycle bface_loop
                
                face_sides = bound(ib)%bfaces(1,j)

                gradqb = zero
                do k = 1,face_sides
                    nk = bound(ib)%bfaces(k + 1,j)
                    gradqb = gradqb + vgradq(:,:,nk)
                end do
                gradqb = gradqb / real(face_sides, p2)

                mu1 = mu(c1)
                mu2 = compute_viscosity(qb(iT))
                muf = half * (mu1 + mu2) ! we do this here because we need it more than once
                if (iflow_type == FLOW_RANS) then
                    trbv1 = turb_var(c1,:)
                    call turb_rhstate(trbv1, ibc_type(ib), trbv2)
                    mutf = calcmut(q1,q2,muf,trbv1,trbv2)
                    ! no elseif needed, we set this to zero before the loop.
                end if

                call visc_flux_boundary(q1,qb,muf,mutf,gradqb,unit_face_normal, &
                                cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                                                              num_flux )

                res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)

#ifdef NANCHECK
                if (any(isnan(res(:,c1)))) then 
                    write (*,*) "nan value present - press [Enter] to continue"
                    read(unit=*,fmt=*)
                end if
#endif

            end do bface_loop

        end do boundary_loop

    end subroutine compute_residual_flow

end module residual