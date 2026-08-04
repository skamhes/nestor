module reorder

    implicit none


    private

    public reorder_rcm

    contains

    subroutine reorder_rcm

        use common , only : p2

        use config , only : rcm_verbosity
        
        use grid , only : cc_data_type, bgrid_type, ncells, cell, face, nfaces, nb, bound, gcell, build_ghost_cells, &
                          face_centroid, face_nrml, face_nrml_mag

        use sort_routines , only : inserstion_sort_ind

        implicit none

        type fc
            integer, dimension(:), allocatable :: iface    ! index of faces
            integer                            :: nf       ! number of faces
        end type fc

        type(cc_data_type), dimension(:)  , pointer :: rcm_cell ! reordered array
        integer,            dimension(:,:), pointer :: rcm_face ! reordered face array
        real(p2),           dimension(:,:), pointer :: rcm_fc, rcm_fn ! other reordered face array
        real(p2),           dimension(:)  , pointer :: rcm_fa

        integer :: min_deg, imd, max_deg ! index of the cell with the minimum degree (# of neighbors)
        integer :: i, j, ib, ihead, iqueue, phead
        integer :: nadj, cn, iadj
        integer :: cL, cR, cLq, cRq
        integer :: fn, nface_loc

        integer, dimension(:,:), allocatable :: adj

        integer, dimension(ncells) :: queue
        
        integer, dimension(ncells) :: c2q ! convert input array to queue
        type(fc),dimension(ncells) :: c2f ! converts cell faces to face array index

        write(*,*) "Reordering mesh using Reverse Cuthill-Mckee"

        c2q = 0

        allocate(rcm_cell(ncells))

        min_deg =  1000000
        max_deg = -1000000

        do i = 1,ncells
            if (cell(i)%nnghbrs < min_deg) then
                min_deg = cell(i)%nnghbrs
                imd     = i
            end if
            max_deg = max(max_deg,cell(i)%nnghbrs)
        end do

        allocate(adj(2,max_deg)) ! adj(:,i) = (/degree, icell/)

        ! starting at a low order cell should yield a better partition.

        ! ihead  : points to the current cell in the old array (will jump around)
        ! iqueue : points to where the current cell is going in the new array (starts at ncells and dec by 1)
        
        ihead  = ncells
        iqueue = ncells

        queue(iqueue) = imd
        c2q(imd)      = iqueue  

        queue_loop : do while(iqueue > 1)
            phead = queue(ihead)

            nadj = cell(phead)%nnghbrs

            iadj = 0

            ! collect adjacent cells that are unset
            do i = 1,nadj
                cn = cell(phead)%nghbr(i)
                if (c2q(cn) == 0) then ! if not set
                    iadj = iadj + 1
                    adj(1,iadj) = cell(cn)%nnghbrs
                    adj(2,iadj) = cn
                end if
            end do
            nadj = iadj

            if (nadj == 0) then
                ihead = ihead - 1
                cycle queue_loop
            elseif(nadj > 1) then 
                ! sort from lowest to highest degree
                call inserstion_sort_ind(nadj, adj)
            end if

            do i = 1,nadj
                iqueue = iqueue - 1
                cn = adj(2,i)
                queue(iqueue) = cn
                c2q(cn) = iqueue
            end do
            ihead = ihead - 1

        end do queue_loop

        deallocate(adj)

        ! Build new cell array

        do i = 1,ncells
            rcm_cell(i) = cell(queue(i))
            do j = 1,rcm_cell(i)%nnghbrs
                cn = rcm_cell(i)%nghbr(j)
                cn = c2q(cn)
                rcm_cell(i)%nghbr(j) = cn
            end do
        end do

        ! This is gonna vomit a bunch of lines to the console.  So only do it if you really mean it.  It's mainly for making nice
        ! pictures from small test grids...
        if (rcm_verbosity == 1001) then
            write(*,*) "Original cell structure:"
            call write_struct(ncells, cell)
            write(*,*) "New cell structure:"
            call write_struct(ncells, rcm_cell)
        endif

        deallocate(cell)

        cell => rcm_cell

        nullify(rcm_cell)

        ! update face array

        allocate(rcm_face(2,ncells)) ! only store the cells.
        
        do i = 1,ncells
            allocate( c2f(i)%iface(cell(i)%nnghbrs) )
            c2f(i)%nf = 0
        end do

        

        do i = 1,nfaces
            ! first make sure the face is pointing to the updated cell indices
            cL  = face(1,i)
            cR  = face(2,i)
            cLq = c2q(cL)
            cRq = c2q(cR)

            face(1,i) = cLq
            face(2,i) = cRq

            ! now grab that data to make the association struct
            cn = minval(face(1:2,i)) ! smaller number face
            c2f(cn)%nf = c2f(cn)%nf + 1
            c2f(cn)%iface(c2f(cn)%nf) = i
        end do

        allocate(rcm_face(2,nfaces))
        allocate(rcm_fc(3,nfaces))
        allocate(rcm_fn(3,nfaces))
        allocate(rcm_fa(nfaces))
        
        nface_loc = 0
        do i = 1,ncells
            do j = 1,c2f(i)%nf
                fn = c2f(i)%iface(j)
                nface_loc = nface_loc + 1
                rcm_face(1:2,nface_loc) = face(1:2,fn)
                rcm_fc(:,nface_loc) = face_centroid(:,fn)
                rcm_fn(:,nface_loc) = face_nrml(:,fn)
                rcm_fa(nface_loc)   = face_nrml_mag(fn)
            end do
            deallocate( c2f(i)%iface )
        end do

        deallocate(face, face_centroid, face_nrml, face_nrml_mag)
        face          => rcm_face
        face_centroid => rcm_fc
        face_nrml     => rcm_fn
        face_nrml_mag => rcm_fa

        nullify(rcm_face, rcm_fc, rcm_fn, rcm_fa)
        
        ! Update boundary array
        do ib = 1,nb
            do i = 1,bound(ib)%nbfaces
                cn = bound(ib)%bcell(i)
                cn = c2q(cn)

                bound(ib)%bcell(i) = cn
            end do
        end do

        if (associated(gcell)) then ! this is a bit sloppy but it works...
            deallocate(gcell)
            call build_ghost_cells
        end if

    end subroutine reorder_rcm

    ! Some diagnostic routines:

    subroutine write_struct(ncells, cell)

        use common , only : p2

        use grid , only : cc_data_type

        implicit none

        integer,                                   intent(in   ) :: ncells
        type(cc_data_type), dimension(:), pointer, intent(inout) :: cell

        integer :: i, j, cn
        character, dimension(ncells + 2) :: line
        integer :: bwidth

        real(p2) :: rbwidth

        bwidth = 0

        do i = 1,ncells
            line(:) = " "
            line(1) = "|"
            line(ncells + 2) = "|"
            do j = 1,cell(i)%nnghbrs
                cn = cell(i)%nghbr(j) + 1
                line(cn) = "x"
            end do
            line(i+1) = "x"
            write(*,*) line
            bwidth = bwidth + max(abs(maxval(cell(i)%nghbr)-i),abs(minval(cell(i)%nghbr)-i))
        end do

        rbwidth = real(bwidth,p2) / real(ncells,p2)

        write(*,*) "Average bandwidth: ", rbwidth
        write(*,*)


    end subroutine write_struct


end module reorder