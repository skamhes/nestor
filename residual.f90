module residual

    implicit none 
    
    public :: compute_residual

    contains 

    subroutine compute_residual

        use common          , only : p2, zero, half, one, two

        use config          , only : method_inv_flux, accuracy_order, use_limiter, &
                                     eps_weiss_smith, method_inv_jac, turbulence_type

        use grid            , only : ncells, cell,  &
                                     nfaces, face,  &
                                     nb,     bound, &
                                     bc_type,       &
                                     face_nrml,     &
                                     face_nrml_mag, &
                                     face_centroid
        
        use solution        , only : res, q, ccgradq, vgradq, wsn, q2u, phi, ur2, compute_uR2

        use interface       , only : interface_flux

        use limiter         , only : compute_limiter

        use bc_states       , only : get_right_state

        use gradient        , only : compute_gradient

        use viscous_flux    , only : visc_flux_boundary, visc_flux_internal

        implicit none

        ! Grid Vars
        real(p2)                    :: xm, ym, zm
        integer                     :: c1, c2,  v1, v2, v3
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2

        ! Flow variables
        real(p2), dimension(5)      :: u1, u2, q1, q2
        real(p2), dimension(3,5)    :: gradq1, gradq2, gradqb
        real(p2), dimension(5)      :: num_flux
        real(p2), dimension(5)      :: qb
        real(p2)                    :: wave_speed
        real(p2)                    :: phi1, phi2
        real(p2)                    :: uR21, uR22

        ! Misc int/counters
        integer                     :: i, os
        integer                     :: j, ib, ix, iu, ii
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
        if (accuracy_order == 2 .OR. trim(turbulence_type) == 'laminar') then
            call compute_gradient(0) ! For now we are just using unweighted gradients
        endif

        if (use_limiter) call compute_limiter

        ! Compute low mach reference velocity
        if(trim(method_inv_flux)=="roe_lm_w" .OR. trim(method_inv_jac)=='roe_lm_w') call compute_uR2


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
            if(trim(method_inv_flux)=="roe_lm_w") then
                uR21 = ur2(c1)
                uR22 = ur2(c2)
            endif
            call interface_flux(          q1,       q2   , & !<- Left/right states
                                      gradq1,      gradq2, & !<- Left/right gradients
                                         unit_face_normal, & !<- unit face normal
                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, & !<- Left  cell centroid
                    cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, & !<- Right cell centroid
                                       face_centroid(1,i), &
                                       face_centroid(2,i), &
                                       face_centroid(3,i), & !<- face midpoint
                                        phi1,        phi2, & !<- Limiter functions
                                               uR21, uR22, &
                                     num_flux, wave_speed ) !<- Output
            ! ur21 & 2 get passed regardless if they have been assigned values.  This is ok since they are only used if they've been
            ! assigned.  Is this sloppy? Maybe?
 
            res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)
            wsn(c1)   = wsn(c1) + wave_speed*face_nrml_mag(i)
            
            res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)
            wsn(c2)   = wsn(c2) + wave_speed*face_nrml_mag(i)

            if ( trim(turbulence_type) == 'inviscid' ) cycle loop_faces

            ! Viscous flux
            call visc_flux_internal(q1,q2,gradq1,gradq2,unit_face_normal,  &
                                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                                    cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, &
                                                                   num_flux)

            res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)

            res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)

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
                call get_right_state(q1, unit_face_normal, bc_type(ib), qb)
                if ( accuracy_order == 2 ) then
                    gradq1 = ccgradq(1:3,1:5,c1)
                else
                    gradq1 = zero
                endif
                
                if(trim(method_inv_flux)=="roe_lm_w") then
                    uR21 = ur2(c1)
                    ! At some point I will probably want to revisit this boundary treatment.  I'm not sure if it will cause issues 
                    ! for now...
                    uR22 = ( min( max( eps_weiss_smith,sqrt(qb(2)**2 + qb(3)**2 + qb(4)**2) ), one) )**2
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
                                               uR21, uR22, &
                                        num_flux, wave_speed  )

                res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)
                wsn(c1)   = wsn(c1) + wave_speed * bound(ib)%bface_nrml_mag(j)

                if ( trim(turbulence_type) == 'inviscid' ) cycle bface_loop
                
                face_sides = bound(ib)%bfaces(1,j)

                gradqb = zero
                do k = 1,face_sides
                    nk = bound(ib)%bfaces(k + 1,j)
                    gradqb = gradqb + vgradq(:,:,nk)
                end do
                gradqb = gradqb / real(face_sides, p2)

                call visc_flux_boundary(q1,qb,gradqb,unit_face_normal, &
                                cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                                                              num_flux )

                res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)

            end do bface_loop

        end do boundary_loop

    end subroutine compute_residual

end module residual