module jacobian

    use common          , only : p2

    use solution        , only : jacobian_type

    implicit none

    public :: compute_jacobian


    contains

    subroutine compute_jacobian

        use common              , only : p2, zero, half, one

        use utils               , only : iturb_type, TURB_INVISCID, ibc_type, ilsq_stencil, LSQ_STENCIL_WVERTEX

        use grid                , only : ncells, nfaces, & 
                                         face, cell, &
                                         face_nrml_mag, face_nrml, &
                                         bound, nb

        use solution            , only : q, gamma, gammamo, gmoinv, dtau, jac, &
                                         kth_nghbr_of_1, kth_nghbr_of_2, ccgradq, vgradq, compute_primative_jacobian

        use interface_jacobian  , only : interface_jac

        use bc_states           , only : get_right_state

        use direct_solve        , only : gewp_solve

        use ad_viscous_flux     , only : visc_flux_boundary_ddt, visc_flux_internal_ddt

        implicit none
        ! Local Vars
        integer                     :: c1, c2, i, k, ib, idestat, j, nk
        real(p2), dimension(3)      :: unit_face_nrml, bface_centroid
        real(p2), dimension(5)      :: qb, q1
        real(p2), dimension(3,5)    :: gradq1, gradq2, gradqb
        real(p2), dimension(5,5)    :: dFnduL, dFnduR
        real(p2)                    :: face_mag

        real(p2), dimension(3,5) :: dummy1, dummy2

        real(p2), dimension(5,5)    :: preconditioner
        
        integer                     :: face_sides

        ! Initialize jacobian terms
        do i = 1,ncells
            jac(i)%diag = zero
            jac(i)%off_diag = zero
            jac(i)%diag_inv = zero
        end do

        ! Loop Faces
        loop_faces : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)

            unit_face_nrml = face_nrml(1:3,i)
            face_mag       = face_nrml_mag(i)

            ! Compute the flux Jacobian for given q1 and q2
            call interface_jac( q(:,c1), q(:,c2), unit_face_nrml, dFnduL, dFnduR)

            ! Add to diagonal term of C1
            jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag
            ! get neighbor index k for cell c1
            k = kth_nghbr_of_1(i)
            ! add to off diagonal neighbor k for cell c1
            jac(c1)%off_diag(:,:,k) = jac(c1)%off_diag(:,:,k) + dFnduR * face_mag

            ! Subtract terms from c2
            jac(c2)%diag            = jac(c2)%diag            - dFnduR * face_mag
            k = kth_nghbr_of_2(i)
            jac(c2)%off_diag(:,:,k) = jac(c2)%off_diag(:,:,k) - dFnduL * face_mag

            if ( iturb_type == TURB_INVISCID ) cycle loop_faces

            gradq1 = ccgradq(1:3,1:5,c1)
            gradq2 = ccgradq(1:3,1:5,c2)

            call visc_flux_internal_ddt(q(:,c1),q(:,c2),gradq1,gradq2,unit_face_nrml, &
                                               cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                                               cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, &
                                                                        dFnduL, dFnduR)
            
            jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag
            ! get neighbor index k for cell c1
            k = kth_nghbr_of_1(i)
            ! add to off diagonal neighbor k for cell c1
            jac(c1)%off_diag(:,:,k) = jac(c1)%off_diag(:,:,k) + dFnduR * face_mag

            ! Subtract terms from c2
            jac(c2)%diag            = jac(c2)%diag            - dFnduR * face_mag
            k = kth_nghbr_of_2(i)
            jac(c2)%off_diag(:,:,k) = jac(c2)%off_diag(:,:,k) - dFnduL * face_mag

        end do loop_faces

        bound_loop : do ib = 1,nb
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(i)
                
                bface_centroid = bound(ib)%bface_center(:,i)
                unit_face_nrml = bound(ib)%bface_nrml(:,i)
                face_mag       = bound(ib)%bface_nrml_mag(i)

                q1 = q(:,c1)
                
                call get_right_state(q1, unit_face_nrml, ibc_type(ib), dummy1, qb, dummy2)

                call interface_jac( q1, qb, unit_face_nrml, dFnduL, dFnduR)
                
                ! We only have a diagonal term to add
                jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag

                if ( iturb_type == TURB_INVISCID ) cycle bfaces_loop

                face_sides = bound(ib)%bfaces(1,i)

                if (ilsq_stencil == LSQ_STENCIL_WVERTEX) then
                    gradqb = zero
                    do k = 1,face_sides
                        nk = bound(ib)%bfaces(k + 1,i)
                        gradqb = gradqb + vgradq(:,:,nk)
                    end do
                    gradqb = gradqb / real(face_sides, p2)
                else ! ilsq_stencil == LSQ_STENCIL_NN
                    gradqb = ccgradq(1:3,1:5,c1)
                endif
                
                call visc_flux_boundary_ddt(q1,qb,gradqb,unit_face_nrml, &
                                  cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                  bface_centroid(1),bface_centroid(2),bface_centroid(3), &
                                                           dFnduL, dFnduR)

                ! We only have a diagonal term to add
                jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag
                
            end do bfaces_loop
        
        end do bound_loop

        ! Now we need to add the pseudo time vol/dtau to the diagonal term along with the jacobian
        ! DQ/DW and generate the inverse diagonal block
        do i = 1,ncells
            preconditioner = compute_primative_jacobian(q(:,i))

            jac(i)%diag = jac(i)%diag + (cell(i)%vol/dtau(i))*preconditioner
            
            ! Invert the diagonal
            idestat = 0
            !                A                 dim  A^{-1}           error check
            call gewp_solve( jac(i)%diag(:,:), 5  , jac(i)%diag_inv, idestat    )
             !  Report errors
            if (idestat/=0) then
                write(*,*) " Error in inverting the diagonal block... Stop"
                write(*,*) "  Cell number = ", i
                do k = 1, 5
                    write(*,'(12(es8.1))') ( jac(i)%diag(k,j), j=1,5 )
                end do
                stop
            endif
        end do
    end subroutine compute_jacobian

end module jacobian

! for later:
