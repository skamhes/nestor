module residual

    implicit none 
    
    public :: compute_residual

    contains 

    subroutine compute_residual

        use common          , only : p2, zero, half, one, two, three_half, three, ix, iy, iz

        use config          , only : method_inv_flux, accuracy_order, use_limiter, mms_include, sutherland_constant, &
                                     reference_temp, M_inf, Re_inf, pr

        use utils           , only : iturb_type, TURB_INVISCID, ilsq_stencil, LSQ_STENCIL_WVERTEX, LSQ_STENCIL_NN

        use grid            , only : ncells, cell,  &
                                     nfaces, face,  &
                                     nb,     bound, &
                                     face_nrml,     &
                                     face_nrml_mag, &
                                     face_centroid

        use utils           , only : ibc_type
        
        use solution        , only : res, q, ccgradq, vgradq, wsn, q2u, phi, iu, iv, iw, iT, gamma, fcgradq

        use interface       , only : interface_flux, reconstruct_flow

        use limiter         , only : compute_limiter

        use bc_states       , only : get_right_state

        use gradient        , only : compute_gradient

        use viscous_flux    , only : visc_flux_boundary, visc_flux_internal

        use mms

        implicit none

        ! Grid Vars
        integer                     :: c1, c2
        real(p2), dimension(3)      :: unit_face_normal, bface_centroid

        ! Flow variables
        real(p2), dimension(5)      :: q1, q2, qL, qR
        real(p2), dimension(3,5)    :: gradq1, gradq2, gradqb
        real(p2), dimension(5)      :: num_flux
        real(p2), dimension(5)      :: qb
        real(p2)                    :: wave_speed
        real(p2)                    :: phi1, phi2
        real(p2)                    :: xc,yc,zc, xc2, yc2,zc2
        real(p2)                    :: fxc,fyc,fzc, dxc2, dyc2,dzc2
        real(p2)                    :: C0, T, mu
        real(p2)                    :: tau_xx, tau_xy, tau_xz
        real(p2)                    :: tau_yx, tau_yy, tau_yz
        real(p2)                    :: tau_zx, tau_zy, tau_zz
        real(p2)                    :: tau_nx, tau_ny, tau_nz
        real(p2)                    :: tau_nV, qx, qy, qz, qn
        real(p2)                    :: ux, uy, uz, vx, vy, vz, wx, wy, wz

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
        if (accuracy_order == 2 .OR. iturb_type > TURB_INVISCID) then
            call compute_gradient(0) ! For now we are just using unweighted gradients
        endif

        ! Only needs to be set once.
        gradq1 = zero
        gradq2 = zero

        phi1 = one
        phi2 = one
        
        if (use_limiter) call compute_limiter

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
                ! Limiters
                if (use_limiter) then
                    phi1 = phi(c1)
                    phi2 = phi(c2)
                endif
                call reconstruct_flow((/cell(c1)%xc, cell(c1)%yc, cell(c1)%zc/)    , &
                    (/face_centroid(1,i), face_centroid(2,i), face_centroid(3,i) /), &
                                                                phi1, q1, gradq1, qL )
                call reconstruct_flow((/cell(c2)%xc, cell(c2)%yc, cell(c2)%zc/)    , &
                    (/face_centroid(1,i), face_centroid(2,i), face_centroid(3,i) /), &
                                                                phi2, q2, gradq2, qR )
            else
                qL = q1
                qR = q2
            endif
            ! Face normal
            unit_face_normal = face_nrml(1:3,i)
            
            call interface_flux(          qL,       qR   , & !<- Left/right states
                                         unit_face_normal, & !<- unit face normal
                                     num_flux, wave_speed  ) !<- Output
            
            if (mms_include(1)) then
                res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)
                wsn(c1)   = wsn(c1) + wave_speed * face_nrml_mag(i)
                
                res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)
                wsn(c2)   = wsn(c2) + wave_speed * face_nrml_mag(i)
            endif

            if ( iturb_type == TURB_INVISCID) cycle loop_faces
            if (mms_include(2)) then
                ! Viscous flux
                fxc = face_centroid(1,i)
                fyc = face_centroid(2,i)
                fzc = face_centroid(3,i)
                call fMMS(fxc, fyc, fzc, qb, qL, gradqb)
                call visc_flux_internal(q1,q2,ccgradq(:,:,c1),ccgradq(:,:,c2), &
                                                             unit_face_normal, &
                                        cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                                        cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, &
                                                        fcgradq(:,:,i),num_flux)
                fcgradq(1:2,:,i) = fcgradq(1:2,:,i) - gradqb(1:2,:)
                ! call visc_flux_internal(q1,q2,gradqb,gradqb, &
                !                                             unit_face_normal, &
                !                         cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                !                         cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, &
                !                                                     num_flux)
                res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)

                res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)
            endif
            if(mms_include(3)) then
                fxc = face_centroid(1,i)
                fyc = face_centroid(2,i)
                fzc = face_centroid(3,i)
                call fMMS(fxc, fyc, fzc, qb, q2, gradqb)
                C0= sutherland_constant/reference_temp
                T = qb(5)
                mu =  M_inf/Re_inf * (one + C0) / (T + C0)*T**(three_half)
                ux = gradqb(ix,iu); vx = gradqb(ix,iv); wx = gradqb(ix,iw);
                uy = gradqb(iy,iu); vy = gradqb(iy,iv); wy = gradqb(iy,iw);
                uz = gradqb(iz,iu); vz = gradqb(iz,iv); wz = gradqb(iz,iw);

                tau_xx = (two/three) * mu * ( two*ux - vy - wz)
                tau_yy = (two/three) * mu * (-ux + two*vy - wz)
                tau_zz = (two/three) * mu * (-ux - vy + two*wz)
                
                tau_xy = mu * (uy + vx); tau_yx = tau_xy;
                tau_xz = mu * (wx + uz); tau_zx = tau_xz;
                tau_yz = mu * (wy + wz); tau_zy = tau_yz;

                tau_nx = -( tau_xx * unit_face_normal(ix) + &
                            tau_xy * unit_face_normal(iy) + &
                            tau_xz * unit_face_normal(iz) )

                tau_ny = -( tau_yx * unit_face_normal(ix) + &
                            tau_yy * unit_face_normal(iy) + &
                            tau_yz * unit_face_normal(iz) )

                tau_nz = -( tau_zx * unit_face_normal(ix) + &
                            tau_zy * unit_face_normal(iy) + &
                            tau_zz * unit_face_normal(iz) )

                tau_nV = tau_nx * qb(iu) + tau_ny * qb(iv) + tau_nz * qb(iw)

                qx     = - mu*gradqb(ix,iT)/(pr*(gamma-one))
                qy     = - mu*gradqb(iy,iT)/(pr*(gamma-one))
                qz     = - mu*gradqb(iz,iT)/(pr*(gamma-one))
                qn     = qx * unit_face_normal(ix) + &
                         qy * unit_face_normal(iy) + &
                         qz * unit_face_normal(iz) 

                num_flux(1) = zero
                num_flux(iu) = tau_nx
                num_flux(iv) = tau_ny
                num_flux(iw) = tau_nz
                
                ! note: + tau_nV because the negative is already baked into the tau_nx/y/z lines
                num_flux(iT) = tau_nV + qn

                res(:,c1) = res(:,c1) + num_flux * face_nrml_mag(i)

                res(:,c2) = res(:,c2) - num_flux * face_nrml_mag(i)
            endif
        end do loop_faces

        !!!!!!!!!!!!!!!!!!!!!!!!V!!!!!!!!!!!!!!!!!!!
        boundary_loop : do ib = 3,nb
            bface_loop : do j = 1,bound(ib)%nbfaces
                bface_centroid = bound(ib)%bface_center(:,j)
                
                c1 = bound(ib)%bcell(j)

                unit_face_normal = bound(ib)%bface_nrml(:,j)

                q1 = q(:,c1)
                
                ! Get the right hand state (weak BC!)
                if ( accuracy_order == 2 ) then
                    gradq1 = ccgradq(1:3,1:5,c1)
                    if (use_limiter) then
                        phi1 = phi(c1)
                        phi2 = phi(c2)
                    endif
                    call reconstruct_flow((/cell(c1)%xc, cell(c1)%yc, cell(c1)%zc/)    , &
                           (/bface_centroid(1), bface_centroid(2), bface_centroid(3) /), &
                                                                    phi1, q1, gradq1, qL )
                else
                    qL = q1
                    qR = q2
                endif

                fxc  = bound(ib)%bface_center(1,j)
                fyc  = bound(ib)%bface_center(2,j)
                fzc  = bound(ib)%bface_center(3,j)
                call get_right_state(qL, (/fxc,fyc,fzc/), unit_face_normal, ibc_type(ib), qb)

                call interface_flux(          qL,      qb, & !<- Left/right states
                                         unit_face_normal, & !<- unit face normal
                                        num_flux, wave_speed  )
                if (mms_include(1)) then
                    res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)
                    wsn(c1)   = wsn(c1) + wave_speed * bound(ib)%bface_nrml_mag(j)
                endif

                if ( iturb_type == TURB_INVISCID ) cycle bface_loop
                if (mms_include(2)) then

                    if (ilsq_stencil == LSQ_STENCIL_WVERTEX) then
                        face_sides = bound(ib)%bfaces(1,j)
                        gradqb = zero
                        do k = 1,face_sides
                            nk = bound(ib)%bfaces(k + 1,j)
                            gradqb = gradqb + vgradq(:,:,nk)
                        end do
                        gradqb = gradqb / real(face_sides, p2)
                    else ! ilsq_stencil == LSQ_STENCIL_NN
                        gradqb = ccgradq(1:3,1:5,c1)
                    endif
                    
                    xc   = cell(c1)%xc
                    yc   = cell(c1)%yc
                    zc   = cell(c1)%zc
                    dxc2 = fxc - xc
                    dyc2 = fyc - yc
                    dzc2 = fzc - zc
                    xc2  = fxc + dxc2
                    yc2  = fyc + dyc2
                    zc2  = fzc + dzc2
                    call get_right_state(q1, (/xc2,yc2,zc2/), unit_face_normal, ibc_type(ib), qb)

                    call fMMS(fxc, fyc, fzc, qL, qR, gradq1)
                    ! q1 = qb
                    call visc_flux_boundary(q1,qb,gradqb,unit_face_normal, &
                                    cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                                                            xc2,yc2,zc2, &
                                                            gradq2, num_flux )
                    
                    fcgradq(1:2,:,1) = fcgradq(1:2,:,1) - gradq1(1:2,:)
                    res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)
                end if
                if(mms_include(3)) then
                    call fMMS(fxc, fyc, fzc, qb, q2, gradqb)
                    C0= sutherland_constant/reference_temp
                    T = qb(5)
                    mu =  M_inf/Re_inf * (one + C0) / (T + C0)*T**(three_half)
                    ux = gradqb(ix,iu); vx = gradqb(ix,iv); wx = gradqb(ix,iw);
                    uy = gradqb(iy,iu); vy = gradqb(iy,iv); wy = gradqb(iy,iw);
                    uz = gradqb(iz,iu); vz = gradqb(iz,iv); wz = gradqb(iz,iw);

                    tau_xx = (two/three) * mu * ( two*ux - vy - wz)
                    tau_yy = (two/three) * mu * (-ux + two*vy - wz)
                    tau_zz = (two/three) * mu * (-ux - vy + two*wz)
                    
                    tau_xy = mu * (uy + vx); tau_yx = tau_xy;
                    tau_xz = mu * (wx + uz); tau_zx = tau_xz;
                    tau_yz = mu * (wy + wz); tau_zy = tau_yz;

                    tau_nx = -( tau_xx * unit_face_normal(ix) + &
                                      tau_xy * unit_face_normal(iy) + &
                                      tau_xz * unit_face_normal(iz) )
                    tau_ny = -( tau_yx * unit_face_normal(ix) + &
                                      tau_yy * unit_face_normal(iy) + &
                                      tau_yz * unit_face_normal(iz) )
                    tau_nz = -( tau_zx * unit_face_normal(ix) + &
                                      tau_zy * unit_face_normal(iy) + &
                                      tau_zz * unit_face_normal(iz) )

                    tau_nV = tau_nx * qb(iu) + tau_ny * qb(iv) + tau_nz * qb(iw)

                    qx     = - mu*gradqb(ix,iT)/(pr*(gamma-one))
                    qy     = - mu*gradqb(iy,iT)/(pr*(gamma-one))
                    qz     = - mu*gradqb(iz,iT)/(pr*(gamma-one))
                    qn     = qx * unit_face_normal(ix) + &
                             qy * unit_face_normal(iy) + &
                             qz * unit_face_normal(iz) 

                    num_flux(1) = zero
                    num_flux(iu) = tau_nx
                    num_flux(iv) = tau_ny
                    num_flux(iw) = tau_nz

                    ! note: + tau_nV because the negative is already baked into the tau_nx/y/z lines
                    num_flux(iT) = tau_nV + qn

                            
                    res(:,c1) = res(:,c1) + num_flux * bound(ib)%bface_nrml_mag(j)
                end if
            end do bface_loop

        end do boundary_loop

    end subroutine compute_residual

end module residual