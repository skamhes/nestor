module grid_statists

    use common , only : p2
    implicit none

    public :: cell_aspect_ratio
    public :: compute_aspect_ratio
    public :: init_ar_array
    public :: grid_spacing
    private :: get_face_nsides
    
    real(p2), dimension(:), pointer :: cell_aspect_ratio
    real(p2), dimension(:), pointer :: grid_spacing ! Nominal delta_x, defined as the distance to the closest cell center.

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

    subroutine compute_grid_spacing
        
        use common  , only : p2, zero, one

        use grid    , only : ncells, cell, nfaces, face

        implicit none

        real(p2) :: xc1, yc1, zc1
        real(p2) :: xc2, yc2, zc2
        real(p2) :: dist
        integer  :: c1, c2
        integer  :: iface

        if (.NOT.associated(grid_spacing)) allocate(grid_spacing(ncells))

        grid_spacing = - one

        do iface = 1,ncells
            c1 = face(1,iface)
            c2 = face(2,iface)
            
            xc1 = cell(c1)%xc
            yc1 = cell(c1)%yc
            zc1 = cell(c1)%zc

            xc2 = cell(c2)%xc
            yc2 = cell(c2)%yc
            zc2 = cell(c2)%zc

            dist = sqrt((xc2-xc1)**2 + (yc2-yc1)**2 + (zc2-zc1)**2)

            ! dist >= 0 will always be true.  So grid_spacing < 0 only occurs if the value is unassigned.
            if ( grid_spacing(c1) < zero ) then
                grid_spacing(c1) = dist
            else
                grid_spacing(c1) = min(grid_spacing(c1),dist)
            endif

            if ( grid_spacing(c2) < zero ) then
                grid_spacing(c2) = dist
            else
                grid_spacing(c2) = min(grid_spacing(c2),dist)
            endif

        enddo

    end subroutine compute_grid_spacing
end module grid_statists