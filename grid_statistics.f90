module grid_statists

    use common , only : p2
    implicit none

    public :: cell_aspect_ratio
    public :: compute_aspect_ratio
    public :: init_ar_array
    private :: get_face_nsides
    
    real(p2), dimension(:), pointer :: cell_aspect_ratio

    contains

    subroutine compute_aspect_ratio

        ! A subroutine to compute the cell aspect ratio using the same formula as the STAR-CCM+ Cell Aspect Ratio field function

        use common , only : p2, zero

        use grid , only : ncells, cell, face, nfaces, & 
                          face_centroid, face_nrml_mag, face_nrml

        implicit none
        integer                     :: icell, iface
        real(p2), parameter         :: NDIM = 3.0_p2 ! Don't feel like dealing with circular depndancies right now...
        integer                     :: c1, c2
        real(p2)                    :: xc1, yc1, zc1, xc2, yc2, zc2
        real(p2)                    :: xf,  yf,  zf
        real(p2), dimension(3)      :: fc1, fc2                  ! vector from face centroid to cell centroid
        real(p2)                    :: cvol, cnfaces
        real(p2), dimension(ncells) :: car_dim1, car_dim2


        car_dim1 = zero
        car_dim2 = zero

        do iface = 1,nfaces
            c1 = face(1,iface)
            c2 = face(2,iface)

            xc1 = cell(c1)%xc
            yc1 = cell(c1)%yc
            zc1 = cell(c1)%zc

            xc2 = cell(c2)%xc
            yc2 = cell(c2)%yc
            zc2 = cell(c2)%zc

            xf  = face_centroid(1,iface)
            yf  = face_centroid(2,iface)
            zf  = face_centroid(3,iface)

            fc1(1) = xf - xc1
            fc1(2) = yf - yc1
            fc1(3) = zf - zc1

            fc2(1) = xf - xc2
            fc2(2) = yf - yc2
            fc2(3) = zf - zc2
            
            car_dim1(c1) = car_dim1(c1) + face_nrml_mag(iface)
            car_dim1(c2) = car_dim1(c2) + face_nrml_mag(iface)

            car_dim2(c1) = car_dim2(c1) + abs(dot_product(face_nrml(:,iface),fc1))
            car_dim2(c2) = car_dim2(c2) + abs(dot_product(face_nrml(:,iface),fc2))
            
        end do

        do icell = 1,ncells
            cvol = cell(ncells)%vol
            cnfaces = get_face_nsides(cell(ncells)%nvtx)

            cell_aspect_ratio(icell) = ndim * cnfaces * max(cvol,zero) / (car_dim1(icell) * car_dim2(icell))
            
        end do

    end subroutine compute_aspect_ratio

    subroutine init_ar_array

        use common , only : p2

        use grid , only : ncells

        implicit none

        if (.NOT.associated(cell_aspect_ratio)) allocate(cell_aspect_ratio(ncells))

    end subroutine init_ar_array

    pure function get_face_nsides(nvtx)

        ! return number of faces as a real(p2) because that's what is needed above
        implicit none

        integer, intent(in) :: nvtx
        real(p2)            :: get_face_nsides

        select case(nvtx)
        case(4) ! tet
            get_face_nsides = 4.0_p2
        case(5) ! pyr
            get_face_nsides = 5.0_p2
        case(6) ! prism
            get_face_nsides = 5.0_p2
        case(8) ! hex
            get_face_nsides = 6.0_p2
        case default ! Error
            get_face_nsides = -1.0_p2
        end select
    end function get_face_nsides
end module grid_statists