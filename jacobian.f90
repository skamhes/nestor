for later:
        ! Create kth_nghbr arrays
allocate(kth_nghbr_of_1(nfaces))
allocate(kth_nghbr_of_2(nfaces))

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