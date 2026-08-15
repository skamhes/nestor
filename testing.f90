! Various test routines 
module test_mod

    use common , only : p2

    contains

    subroutine tri_diag

        use linear_solver , only : thomas_sweep_block

        use direct_solve , only  : gewp_solve

        implicit none

        integer,  dimension(:), allocatable     :: R,C
        integer,  dimension(3)                  :: Rline
        real(p2), dimension(:,:,:), allocatable :: V, dinv
        real(p2), dimension(:,:),   allocatable :: b, x, x_solve, lres

        integer :: nq, size, nnz
        integer, dimension(:), allocatable :: iline
        integer :: cj, ci
        integer :: i, j, k

        nq = 1; size = 5; nnz = 4 + 3*(size-2) + size
        allocate(V(nq,nq,nnz), b(nq,size), x(nq,size), x_solve(nq,size))
        allocate(R(2*size + 1) , C(nnz))
        allocate(dinv(nq,nq,size))
        Rline = (/ 1, 1+size, 6+size/)
        V(:,:,1:5)   = reshape((/1._p2, 0.5_p2, .7_p2, 0.1_p2, 0.8_p2/),(/nq,nq,5/))
        V(:,:,6:7)   = reshape(  (/5._p2,1.0_p2/),(/nq,nq,2/))
        V(:,:,8:10)  = reshape(  (/1._p2,4._p2,1._p2/),(/nq,nq,3/))
        V(:,:,11:13) = reshape(  (/1._p2,7._p2,2._p2/),(/nq,nq,3/))
        V(:,:,14:16) = reshape( (/2._p2,4._p2,1._p2/),(/nq,nq,3/))
        V(:,:,17:18) = reshape((/1._p2,6._p2/),(/nq,nq,2/))

        C(1:5) = (/4,5,1,1,2/)
        C(6:7) = (/1,2/)
        R(1:size+1) = (/1,2,3,4,5,6/)
        R(size+2) = 8
        call gewp_solve(V(:,:,size + 1), nq, dinv(:,:,1), k)
        do i = size+2,2*size-1
            R(i+1) = R(i) + 3
            C(R(i):R(i)+2) = i-size + (/-1, 0, 1/)
            call gewp_solve(V(:,:,R(i)+1), nq, dinv(:,:,i-size), k)
        end do
        call gewp_solve(V(:,:,nnz), nq, dinv(:,:,size), k)
        C(R(2*size):R(2*size)+1) = size - (/1, 0/)
        R(2*size+1) = nnz+1

        x(:,:) = reshape((/1.4_p2,2.0_p2,1.5_p2,1.5_p2,2._p2/),(/nq,size/))

        b = 0._p2

        do i = 1,2*size
            ci = mod(i-1,5)+1
            do j = R(i),R(i+1)-1
                cj = C(j)
                b(:,ci) = b(:,ci) + matmul(V(:,:,j),x(:,cj))
            end do
        end do

        allocate(iline(size))
        iline = (/1,2,3,4,5/)
        x_solve = 0._p2
        do i = 1,20
            call thomas_sweep_block(size, iline, nq, V, C, R, Rline, Dinv, -b, x_solve,lres, k)
            ! write(*,*) x
            ! write(*,*) x_solve
        end do

        deallocate(V,b,x,x_solve)
        deallocate(R,C,dinv)
        deallocate(iline)

        ! Test a hand built 2x2 block tridiag with strong diagonal dominance
        nq = 2; size = 5; nnz = 4 + 3*(size-2) + size
        allocate(V(nq,nq,nnz), b(nq,size), x(nq,size), x_solve(nq,size))
        allocate(R(2*size + 1) , C(nnz))
        allocate(dinv(nq,nq,size))
        Rline = (/ 1, 1+size, 6+size/)
        V(:,:,1:3) = reshape(  (/0.4_p2, 2._p2, 0.2_p2,0.8_p2, 2._p2,0.8_p2,1._p2,1._p2, 1._p2,2._p2, 0.1_p2, 2._p2/), (/nq,nq,3/) )
        V(:,:,4:5) = reshape(  (/0.4_p2, 2._p2, 0.2_p2,0.1_p2, 2._p2,0.8_p2,1._p2,1._p2/), (/nq,nq,2/) )
        V(:,:,6:7) = reshape(  (/5._p2,1.0_p2,2._p2,6._p2, 1._p2,0.5_p2,1._p2,2._p2/),(/2,2,2/))
        V(:,:,8:10) = reshape(  (/1._p2,2._p2,1._p2,0.5_p2, 8._p2,0.5_p2,1._p2,5._p2, 0.3_p2, 1._p2, 0.6_p2, 1._p2/),(/2,2,3/))
        V(:,:,11:13) = reshape(  (/0.5_p2,1._p2,2._p2,1._p2, 7._p2,0.2_p2,2._p2,6._p2, 0.4_p2, 2._p2, 0.2_p2, 2._p2/),(/2,2,3/))
        V(:,:,14:16) = reshape( (/1._p2,0.7_p2,1._p2,2._p2, 4._p2,0.8_p2,2._p2,6._p2, 0.5_p2, 2._p2, 0.1_p2, 2._p2/),(/2,2,3/))
        V(:,:,17:18) = reshape((/2._p2,0.8_p2,1._p2,1._p2, 6._p2,0.1_p2,1._p2,4._p2/),(/2,2,2/))

        C(1:5) = (/4,5,1,1,2/)
        C(6:7) = (/1,2/)
        R(1:size+1) = (/1,2,3,4,5,6/)
        R(size+2) = 8
        call gewp_solve(V(:,:,size + 1), nq, dinv(:,:,1), k)
        do i = size+2,2*size-1
            R(i+1) = R(i) + 3
            C(R(i):R(i)+2) = i-size + (/-1, 0, 1/)
            call gewp_solve(V(:,:,R(i)+1), nq, dinv(:,:,i-size), k)
        end do
        call gewp_solve(V(:,:,nnz), nq, dinv(:,:,size), k)
        C(R(2*size):R(2*size)+1) = size - (/1, 0/)
        R(2*size+1) = nnz+1

        x(:,:) = reshape((/1.4_p2,2.0_p2, 1.5_p2,1.5_p2, 2._p2, 0.5_p2, 3._p2,2._p2, 1.3_p2, 1.0_p2/),(/2,5/))

        b = 0._p2

        do i = 1,2*size
            ci = mod(i-1,5)+1
            do j = R(i),R(i+1)-1
                cj = C(j)
                b(:,ci) = b(:,ci) + matmul(V(:,:,j),x(:,cj))
            end do
        end do

        allocate(iline(size))
        iline = (/1,2,3,4,5/)
                
        do i = 1,20
            write(*,*) "iteration:", i
            call thomas_sweep_block(size, iline, nq, V, C, R, Rline, Dinv, -b, x_solve, lres, k)
            write(*,*) "Error: ", sum(abs(x_solve - x))
            ! write(*,*) x_solve
        end do

        write(*,*)
        write(*,*)
        do i = 1,size
            write(*,'(i2)') i
            do j = 1,nq
                write(*,*) x(j,i), x_solve(j,i)
            end do
            write(*,*)
        end do

        deallocate(V,b,x,x_solve)
        deallocate(R,C,dinv)
        deallocate(iline)



    end subroutine tri_diag

    subroutine stri_diag

        use linear_solver , only : thomas_sweep_scalar

        implicit none

        integer,  dimension(:), allocatable     :: R,C
        integer,  dimension(3)                  :: Rline
        real(p2), dimension(:), allocatable :: V, dinv, deltai
        real(p2), dimension(:),   allocatable :: b, x, x_solve, res

        integer :: nq, size, nnz
        integer, dimension(:), allocatable :: iline
        integer :: cj, ci
        integer :: i, j, k, iter
        real(p2) :: l, d, u, um1, lres

        nq = 1; size = 5; nnz = 4 + 3*(size-2) + 5
        allocate(V(nnz), b(size), x(size), x_solve(size), deltai(size))
        allocate(R(2*size + 1) , C(nnz))
        allocate(dinv(size), res(size))
        Rline = (/ 1, 1+size, 6+size/)
        V(1:5)     = (/1._p2, 0.5_p2, .7_p2, 0.1_p2, 0.8_p2/)
        V(5+1:5+2) = (/5._p2,1.0_p2/)
        V(5+3:5+5) = (/1._p2,4._p2,1._p2/)
        V(5+6:5+8) = (/1._p2,7._p2,2._p2/)
        V(5+9:5+11) = (/2._p2,4._p2,1._p2/)
        V(5+12:5+13) = (/1._p2,6._p2/)

        C(1:5) = (/4,5,1,1,2/)
        C(6:7) = (/1,2/)
        R(1:size+1) = (/1,2,3,4,5,6/)
        R(size+2) = 8
        dinv(1) = 1._p2 / V( R(Rline(2)) )
        do i = size+2,2*size-1
            R(i+1) = R(i) + 3
            C(R(i):R(i)+2) = i-size + (/-1, 0, 1/)
            dinv(i-size) = 1._p2 / V(R(i)+1)
        end do
        dinv(size) = 1._p2 / V(nnz)
        C(R(2*size):R(2*size)+1) = size - (/1, 0/)
        R(2*size+1) = nnz+1

        x(:) = (/1.4_p2,2.0_p2,1.5_p2,1.5_p2,2._p2/)

        b = 0._p2

        do i = 1,2*size
            ci = mod(i-1,5)+1
            do j = R(i),R(i+1)-1
                cj = C(j)
                b(ci) = b(ci) + V(j)*x(cj)
            end do
        end do
        
        allocate(iline(size))
        iline = (/1,2,3,4,5/)

        res = b

        x_solve = 0._p2
        
        do iter = 1,2
            write(*,*) "iteration:", iter

            do i = 1,size
                res(i) = b(i)
                do j = R(i),R(i+1)-1
                    cj = C(j)
                    res(i) = res(i) - x_solve(cj) * V(j)
                end do
            end do
            deltai(1) = dinv(1)
            res(1)    = res(1) * deltai(1)
            deltai(1) = V( R(Rline(2))+1 ) * deltai(1)
            write(*,*) "del_",1,"= ", deltai(1), "res=",res(1)

            do i = 2,size-1
                l   = V(R(i+size)  )
                d   = V(R(i+size)+1)
                u   = V(R(i+size)+2)

                deltai(i) = 1._p2 / (d - l*deltai(i-1))
                res(i)    = (res(i) - l * res(i-1)) * deltai(i)
                deltai(i) = u * deltai(i)
                write(*,*) "del_",i,"= ", deltai(i), "res=",res(i)

            end do

            l   = V(R(2*size)  )
            d   = V(R(2*size)+1)

            deltai(size) = 1._p2 / (d - l*deltai(size-1))
            res(size)    = (res(size) - l * res(size-1)) * deltai(size)
            write(*,*) "del_",i,"= ", deltai(i), "res=", res(size)

            x_solve(size) = res(size)
            write(*,*) "x_solve_",size,"=",x_solve(size), " x_",size,"=",x(size)
            do i = size-1,1,-1
                x_solve(i) = res(i) - deltai(i)*x_solve(i+1)
                write(*,*) "x_solve_",i,"=",x_solve(i), " x_",i,"=",x(i)
            end do
            write(*,*) "Error: ", sum(abs(x-x_solve))
        end do

        x_solve = 0._p2

        do iter = 1,20
            write(*,*) "iteration:", iter
            lres = 0._p2
            call thomas_sweep_scalar(size, iline, V, C, R, Rline, Dinv, -b, x_solve,lres, k)
            write(*,*) lres
        end do
        
        stop
    end subroutine stri_diag
endmodule test_mod




program testing

    use test_mod

    call stri_diag

    call tri_diag
end program testing