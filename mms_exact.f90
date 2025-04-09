module mms_exact

    use common , only : p2

    implicit none

    contains

    subroutine manual_test_mms(cellin,icell)

        use common , only : ix, iy, iz

        use grid , only : nfaces, cc_data_type, face, face_nrml, face_nrml_mag, face_centroid, cell, bound, nb

        use interface , only : interface_flux

        use viscous_flux , only : compute_visc_num_flux, visc_flux_internal, visc_flux_boundary

        use mms

        use solution , only : ccgradq, q, inv_ncells, vgradq

        implicit none 

        type(cc_data_type), intent(in) :: cellin
        integer,            intent(in) :: icell

        integer :: iface, c1, c2, ib, bcell, nk, k
        real(p2) :: fnrmlM
        real(p2), dimension(3) :: n12
        real(p2) :: fx, fy, fz, cx, cy, cz, dummy1
        real(p2) :: xc1, xc2, yc1, yc2, zc1, zc2

        real(p2), dimension(5) :: S_exact, Qcc, Qfc, num_fluxe, resid, num_flux, reside
        real(p2), dimension(5) :: q1, q2
        real(p2), dimension(3,5) :: gradF, dummy35
        real(p2), dimension(3,5) :: gradq1, gradq2, gradqb
        integer :: face_sides

        real(p2) :: txx_xerr, nf_err

        real(p2) :: xc,yc,zc,fxc,fyc,fzc,dxc2,dyc2,dzc2


        dummy35 = 0.0_p2
        resid = 0.0_p2
        reside = resid
        nf_err = 0.0_p2

        write(*,*)
        write(*,*)
        write(*,*)

        floop : do iface = 1,nfaces
            c1 = face(1,iface)
            c2 = face(2,iface)
            if (c1 == icell) then
                fnrmlM = face_nrml_mag(iface)
                n12    =  face_nrml(:,iface)
                fx     = face_centroid(ix,iface)
                fy     = face_centroid(iy,iface)
                fz     = face_centroid(iz,iface)
            elseif (c2 == icell) then
                fnrmlM = face_nrml_mag(iface)
                n12    =  face_nrml(:,iface) ! always pointing out of the cell
                fx     = face_centroid(ix,iface)
                fy     = face_centroid(iy,iface)
                fz     = face_centroid(iz,iface)
            else
                cycle floop
            endif   

            cx = cellin%xc
            cy = cellin%yc
            cz = cellin%zc
            
            call fMMS(cx, cy, cz, Qcc, S_exact)
            call fMMS(fx, fy, fz, Qfc, S_exact, gradF)
            q1 = q(:,c1)
            q2 = q(:,c2)
            gradq1 = ccgradq(:,:,c1)
            gradq2 = ccgradq(:,:,c2)
            xc1 = cell(c1)%xc
            yc1 = cell(c1)%yc
            zc1 = cell(c1)%zc
            xc2 = cell(c2)%xc
            yc2 = cell(c2)%yc
            zc2 = cell(c2)%zc
            call interface_flux(q1,q2,gradq1,gradq2, n12,xc1,yc1,zc1,xc2,yc2,zc2,fx,fy,fz, 1.0_p2, 1.0_p2, &
                                num_flux, dummy1)
            if (icell == c2) num_flux = - num_flux
            
            ! resid = resid + num_flux * fnrmlM
            

            ! exact values
            call compute_visc_num_flux(Qfc,Qfc,gradF,n12,num_fluxe)
            ! write(*,'(a40,3es13.5)') " exact face grad", gradF(:,2)

            ! numerical values
            call visc_flux_internal(q1,q2,gradq1,gradq2,n12,xc1,yc1,zc1,xc2,yc2,zc2, num_flux)
            if (icell == c2) num_flux = - num_flux
            resid = resid + num_flux * fnrmlM
            reside = reside + num_fluxe * fnrmlM
            
            ! nf_err = nf_err + abs(num_flux(2) - num_fluxe(2))
        !     write(*,'(a,5es13.5,a,i5)') " nfe: ", num_fluxe(:), " face: ", iface
        !     write(*,'(a,5es13.5,a,i5)') " nfn: ", num_flux( :), " face: ", iface
        !     write(*,*)
        end do floop

        ! bface
        bloop : do ib = 3,nb
            floop2 : do iface = 1,bound(ib)%nbfaces
                bcell = bound(ib)%bcell(iface)
                if (bcell /= icell) cycle floop2
                fnrmlM = bound(ib)%bface_nrml_mag(iface)
                n12    = bound(ib)%bface_nrml(:,iface)
                fx     = bound(ib)%bface_center(ix,iface)
                fy     = bound(ib)%bface_center(iy,iface)
                fz     = bound(ib)%bface_center(iz,iface)

                face_sides = bound(ib)%bfaces(1,iface)

                q1 = q(:,icell)
                
                ! Get the right hand state (weak BC!)
                xc   = cell(icell)%xc
                yc   = cell(icell)%yc
                zc   = cell(icell)%zc
                fxc  = bound(ib)%bface_center(1,iface)
                fyc  = bound(ib)%bface_center(2,iface)
                fzc  = bound(ib)%bface_center(3,iface)
                dxc2 = fxc - xc
                dyc2 = fyc - yc
                dzc2 = fzc - zc
                xc2  = fxc + dxc2
                yc2  = fyc + dyc2
                zc2  = fzc + dzc2

                gradqb = 0.0_p2
                do k = 1,face_sides
                    nk = bound(ib)%bfaces(k + 1,iface)
                    gradqb = gradqb + vgradq(:,:,nk)
                end do
                gradqb = gradqb / real(face_sides, p2)

                call fMMS(fx, fy, fz, Qfc, S_exact, gradF)

                call compute_visc_num_flux(Qfc,Qfc,gradF,n12,num_fluxe)
                
                call fMMS(xc2, yc2, zc2, q2, S_exact, gradF)

                call visc_flux_boundary(q1,q2,gradqb,n12, &
                       cell(icell)%xc, cell(icell)%yc, cell(icell)%zc, &
                                                          xc2,yc2,zc2, &
                                                              num_flux )

                resid = resid + num_flux * fnrmlM
                reside = reside + num_fluxe * fnrmlM

                ! write(*,'(a,5es13.5,a,2i5)') " nfe: ", num_fluxe(:), " bound, face: ", ib, iface
                ! write(*,'(a,5es13.5,a,2i5)') " nfn: ", num_flux( :), " bound, face: ", ib, iface
                ! write(*,*)
            end do floop2
        end do bloop
        ! write(*,*) gradF(:,1)
        ! write(*,*) gradF(:,2)
        ! write(*,*) gradF(:,3)
        ! write(*,*) gradF(:,4)
        ! write(*,*) gradF(:,5)
        ! write(*,*) num_fluxe




        call fMMS(cx, cy, cz, Qcc, S_exact)
        ! resid = abs(resid - S_exact * cellin%vol)
        txx_xerr = abs(resid(2) - txx_x * cellin%vol) / cellin%vol * inv_ncells

        ! resid = (resid - S_exact * cellin%vol) * inv_ncells / cellin%vol

        ! write(*,*) " -------------- Truncation error (Manual) ------------------"
        ! write(*,'(a25,i5,a,5es13.5)') " rnumeric_L1(TE) atuo @ icell =  ",icell, ":", abs(resid ) / cellin%vol * inv_ncells
        ! write(*,'(a25,i5,a,5es13.5)') "   rexact_L1(TE) atuo @ icell =  ",icell, ":", abs(reside) / cellin%vol * inv_ncells
        ! write(*,'(a,1es13.5)',advance='no') " numeric: txx_x(TE), tyy_y(TE) ", txx_xerr
        ! txx_xerr = abs(resid(3) - tyy_y * cellin%vol) / cellin%vol * inv_ncells
        ! write(*,'(1es13.5)')  txx_xerr

        ! txx_xerr = abs(reside(2) - txx_x * cellin%vol) / cellin%vol * inv_ncells
        ! write(*,'(a,1es13.5)',advance='no') "  exact: txx_x(TE), tyy_y(TE) ", txx_xerr
        ! txx_xerr = abs(reside(3) - tyy_y * cellin%vol) / cellin%vol * inv_ncells
        ! write(*,'(1es13.5)')  txx_xerr
        ! write(*,'(a,1es13.5)') " NF vs NFe(TE): ", nf_err
    end subroutine manual_test_mms
end module mms_exact