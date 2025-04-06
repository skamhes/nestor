module mms_exact

    use common , only : p2

    implicit none

    contains

    subroutine manual_test_mms(cell,icell)

        use common , only : ix, iy, iz

        use grid , only : nfaces, cc_data_type, face, face_nrml, face_nrml_mag, face_centroid

        use interface , only : interface_flux

        use viscous_flux , only : compute_visc_num_flux

        use mms

        implicit none 

        type(cc_data_type), intent(in) :: cell
        integer,            intent(in) :: icell

        integer :: iface, c1, c2
        real(p2) :: fnrmlM
        real(p2), dimension(3) :: n12
        real(p2) :: fx, fy, fz, cx, cy, cz, dummy1

        real(p2), dimension(5) :: S_exact, Qcc, Qfc, num_flux, resid
        real(p2), dimension(3,5) :: gradF, dummy35
        dummy35 = 0.0_p2
        resid = 0.0_p2

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
                n12    = -face_nrml(:,iface) ! always pointing out of the cell
                fx     = face_centroid(ix,iface)
                fy     = face_centroid(iy,iface)
                fz     = face_centroid(iz,iface)
            else
                cycle floop
            endif   

            cx = cell%xc
            cy = cell%yc
            cz = cell%zc
            
            call fMMS(cx, cy, cz, Qcc, S_exact)
            call fMMS(fx, fy, fz, Qfc, S_exact, gradF)

            call interface_flux(Qfc,Qfc,dummy35,dummy35, n12,cx,cy,cz,fx,fy,fz,fx,fy,fz, 1.0_p2, 1.0_p2, &
                                num_flux, dummy1)

            resid = resid + num_flux * fnrmlM

            call compute_visc_num_flux(Qfc,Qfc,gradF,n12,num_flux)

            resid = resid + num_flux * fnrmlM

        end do floop

        resid = abs(resid - S_exact * cell%vol)

        write(*,*) " -------------- Truncation error (Exact face values) ------------------"
        write(*,'(a,5es13.5)') " rE_L1(TE) ", resid
    end subroutine manual_test_mms
end module mms_exact