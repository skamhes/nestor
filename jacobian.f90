module jacobian

    use common          , only : p2

    implicit none

    public :: init_jacobian
    public :: compute_jacobian

    public :: kth_nghbr_of_1, kth_nghbr_of_2
    integer, dimension(:), allocatable :: kth_nghbr_of_1
    integer, dimension(:), allocatable :: kth_nghbr_of_2

    ! Jacobian type and var
    type jacobian_type
        real(p2), dimension(5,5)                :: diag     ! diagonal blocks of Jacobian matrix
        real(p2), dimension(:,:,:), allocatable :: off_diag ! off-diagonal blocks
        real(p2), dimension(5,5)                :: diag_inv ! inverse of diagonal blocks
        real(p2), dimension(5)                  :: RHS      ! Right hand side (b) of the linear system
    end type jacobian_type

    public :: jac
    type(jacobian_type), dimension(:), allocatable :: jac ! jacobian array
    contains

    
    subroutine init_jacobian

        use common          , only : p2

        use grid            , only : nfaces, face, cell, ncells

        use solution        , only : nq

        implicit none

        integer :: i, k
        integer :: c1, c2

        ! Create kth_nghbr arrays
        allocate(kth_nghbr_of_1(nfaces))
        allocate(kth_nghbr_of_2(nfaces))
        allocate(jac           (ncells))

        ! Define kth neighbor arrays
        face_nghbr_loop : do i = 1,nfaces
            c1 = face(1,i)
            c2 = face(2,i)
            ! loop over c1 neighbors to find c2
            do k = 1,cell(c1)%nnghbrs
                if ( c2 == cell(c1)%nghbr(k)) then
                    kth_nghbr_of_1(i) = k ! c2 is the kth neighbor of c1
                end if
            end do
            ! repeat for cell 2
            do k = 1,cell(c2)%nnghbrs
                if ( c1 == cell(c2)%nghbr(k)) then
                    kth_nghbr_of_2(i) = k
                end if
            end do
        end do face_nghbr_loop

        ! allocate jacobian off diagonal arrays
        do i = 1,ncells
            allocate(  jac(i)%off_diag(nq,nq,cell(i)%nnghbrs))
        end do

    end subroutine init_jacobian

    subroutine compute_jacobian

        use common          , only : p2, zero

        use grid            , only : ncells, nfaces, & 
                                     face, cell, &
                                     face_nrml_mag, face_nrml, &
                                     bound, nb, bc_type

        use solution        , only : q

        implicit none
        ! Local Vars
        integer                     :: c1, c2, i, k, ib, idestat, j, os, ii,jj
        real(p2), dimension(3)      :: unit_face_nrml, bface_centroid, ds, d_Cb, ejk
        real(p2), dimension(5)      :: u1, u2, ub, wb, qb, q1
        real(p2), dimension(3,5)    :: gradw1, gradw2, gradwb
        real(p2), dimension(5,5)    :: dFnduL, dFnduR
        real(p2)                    :: face_mag, mag_ds, mag_ejk
        real(p2)                    :: xc1,xc2,yc1,yc2,zc1,zc2

        real(p2), dimension(5,5)    :: preconditioner, pre_inv
        real(p2), dimension(5,5)    :: duLdqL, duRdqR
        real(p2)                    :: theta
        real(p2)                    :: rho_p, rho_T, rho
        real(p2)                    :: H, alpha, beta, lambda, absu

        ! Initialize jacobian terms
        do i = 1,ncells
            jac(i)%diag = zero
            jac(i)%off_diag = zero
            jac(i)%diag_inv = zero
        end do

        ! Loop Faces
        do i = 1,nfaces
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

        end do

        bound_loop : do ib = 1,nb
            bfaces_loop : do i = 1,bound(ib)%nbfaces
                c1 = bound(ib)%bcell(i)
                
                bface_centroid = bound(ib)%bface_center(:,i)
                unit_face_nrml = bound(ib)%bface_nrml(:,i)
                face_mag       = bound(ib)%bface_nrml_mag(i)

                q1 = q(:,c1)
                
                call get_right_state(q1, unit_face_nrml, bc_type(ib), qb)

                call interface_jac( q1, qb, unit_face_nrml, dFnduL, dFnduR)
                
                ! We only have a diagonal term to add
                jac(c1)%diag            = jac(c1)%diag            + dFnduL * face_mag

            end do bfaces_loop
        
        end do bound_loop

                
                
        
    end subroutine compute_jacobian

end module jacobian

! for later:
