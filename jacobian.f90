module jacobian

    use common          , only : p2

    implicit none

    public :: compute_jacobian


    contains

    subroutine compute_jacobian

        use common              , only : p2, zero, half

        use utils               , only : iflow_type, FLOW_INVISCID, FLOW_RANS, ibc_type, ilsq_stencil, LSQ_STENCIL_WVERTEX

        use grid                , only : ncells, nfaces, & 
                                         face, cell, gcell, &
                                         face_nrml_mag, face_nrml, &
                                         bound, nb, gcell

        use solution_vars       , only : q, dtau, jac, kth_nghbr_of_1, kth_nghbr_of_2, ccgradq, vgradq, iT, kth_of_cell, diag_inv

        use solution            , only : compute_primative_jacobian

        use interface_jacobian  , only : interface_jac

        use bc_states           , only : get_right_state

        use direct_solve        , only : gewp_solve

        use ad_viscous_flux     , only : visc_flux_boundary_ddt, visc_flux_internal_ddt

        use turb                , only : turb_var, nturb, calcmut

        use turb_bc             , only : turb_rhstate

        use viscosity           , only : compute_viscosity

        implicit none
        ! Local Vars
        integer                     :: c1, c2, i, ib, idestat, j, nk, k
        integer                     :: ic1, ic2, k1, k2
        real(p2), dimension(3)      :: unit_face_nrml, bface_centroid
        real(p2), dimension(5)      :: qb, q1
        real(p2), dimension(3,5)    :: gradq1, gradq2, gradqb
        real(p2), dimension(5,5)    :: dFnduL, dFnduR
        real(p2)                    :: face_mag
        real(p2)                    :: xc2, yc2, zc2

        real(p2), dimension(3,5)    :: dummy1, dummy2
        real(p2)                    :: mu1, mu2, muf
        real(p2)                    :: mutf
        real(p2), dimension(nturb)  :: trbv1,trbv2

        real(p2), dimension(5,5)    :: preconditioner
        
        integer                     :: face_sides

        ! Initialize jacobian terms
        jac(:,:,:) = zero

        mutf = zero
        trbv1 = zero
        trbv2 = zero
        ! Loop Faces
        loop_faces : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)

            unit_face_nrml = face_nrml(1:3,i)
            face_mag       = face_nrml_mag(i)

            ! Compute the flux Jacobian for given q1 and q2
            call interface_jac( q(:,c1), q(:,c2), unit_face_nrml, dFnduL, dFnduR)

            ! Add to diagonal term of C1
            ic1 = kth_of_cell(c1)
            jac(:,:,ic1) = jac(:,:,ic1) + dFnduL * face_mag
            ! get neighbor index k for cell c1
            k1 = kth_nghbr_of_1(i)
            ! add to off diagonal neighbor k for cell c1
            jac(:,:,k1)   = jac(:,:,k1)   + dFnduR * face_mag

            ! Subtract terms from c2
            ic2 = kth_of_cell(c2)
            jac(:,:,ic2) = jac(:,:,ic2) - dFnduR * face_mag
            k2 = kth_nghbr_of_2(i)
            jac(:,:,k2)   = jac(:,:,k2)   - dFnduL * face_mag

            if ( iflow_type == FLOW_INVISCID ) cycle loop_faces

            gradq1 = ccgradq(1:3,1:5,c1)
            gradq2 = ccgradq(1:3,1:5,c2)

            if (iflow_type == FLOW_RANS) then
                trbv1 = turb_var(c1,:)
                trbv2 = turb_var(c2,:)
            end if

            call visc_flux_internal_ddt(q(:,c1),q(:,c2),gradq1,gradq2,trbv1,trbv2, &
                                                                   unit_face_nrml, &
                                            cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                                            cell(c2)%xc, cell(c2)%yc, cell(c2)%zc, &
                                                                     dFnduL, dFnduR)
            
            ! Add to diagonal term of C1
            ! ic1 = kth_of_cell(c1)
            jac(:,:,ic1) = jac(:,:,ic1) + dFnduL * face_mag
            ! get neighbor index k for cell c1
            ! k = kth_nghbr_of_1(i)
            ! add to off diagonal neighbor k for cell c1
            jac(:,:,k1)   = jac(:,:,k1)   + dFnduR * face_mag

            ! Subtract terms from c2
            ! ic2 = kth_of_cell(c2)
            jac(:,:,ic2) = jac(:,:,ic2) - dFnduR * face_mag
            ! k = kth_nghbr_of_2(i)
            jac(:,:,k2)   = jac(:,:,k2)   - dFnduL * face_mag

        end do loop_faces

        bound_loop : do ib = 1,nb
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(i)
                
                bface_centroid = bound(ib)%bface_center(:,i)
                unit_face_nrml = bound(ib)%bface_nrml(:,i)
                face_mag       = bound(ib)%bface_nrml_mag(i)

                q1 = q(:,c1)

                xc2  = gcell(ib)%xc(i)
                yc2  = gcell(ib)%yc(i)
                zc2  = gcell(ib)%zc(i)
                
                call get_right_state(q1, unit_face_nrml, ibc_type(ib), qb)

                call interface_jac( q1, qb, unit_face_nrml, dFnduL, dFnduR)
                
                ! We only have a diagonal term to add
                ic1 = kth_of_cell(c1)
                jac(:,:,ic1) = jac(:,:,ic1) + dFnduL * face_mag

                if ( iflow_type == FLOW_INVISCID ) cycle bfaces_loop

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

                xc2  = gcell(ib)%xc(i)
                yc2  = gcell(ib)%yc(i)
                zc2  = gcell(ib)%zc(i)
                
                if (iflow_type == FLOW_RANS) then
                    trbv1 = turb_var(c1,:)
                    call turb_rhstate(trbv1, ibc_type(ib), trbv2)
                end if
                call visc_flux_boundary_ddt(q1,qb,gradqb,trbv1,trbv2, &
                                                      unit_face_nrml, &
                               cell(c1)%xc, cell(c1)%yc, cell(c1)%zc, &
                                                         xc2,yc2,zc2, &
                                                        dFnduL, dFnduR)

                ! We only have a diagonal term to add
                ! ic1 = kth_of_cell(c1)
                jac(:,:,ic1) = jac(:,:,ic1) + dFnduL * face_mag
                
            end do bfaces_loop
        
        end do bound_loop

        ! Now we need to add the pseudo time vol/dtau to the diagonal term along with the jacobian
        ! DQ/DW and generate the inverse diagonal block
        do i = 1,ncells
            preconditioner = compute_primative_jacobian(q(:,i))

            ic1 = kth_of_cell(c1)
            jac(:,:,ic1) = jac(:,:,ic1) + (cell(i)%vol/dtau(i))*preconditioner
            
            ! Invert the diagonal
            idestat = 0
            !                A                 dim  A^{-1}           error check
            call gewp_solve( jac(:,:,ic1), 5  , diag_inv(:,:,i), idestat    )
             !  Report errors
            if (idestat/=0) then
                write(*,*) " Error in inverting the diagonal block... Stop"
                write(*,*) "  Cell number = ", i
                do k = 1, 5
                    write(*,'(12(es8.1))') ( jac(k,j,ic1), j=1,5 )
                end do
                stop
            endif
        end do
    end subroutine compute_jacobian

end module jacobian

! for later:
