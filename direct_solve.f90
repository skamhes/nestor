module direct_solve

    use common , only : p2, zero, one, two

    implicit none

    public :: qr_factorization

    contains

    !****************************************************************************
    ! ------------------ QR Factorization ---------------------
    !
    !  This subroutine solves the LSQ problem: A*x=b, A=mxn matrix.
    !
    !  IN :       a = (m x n) LSQ matrix.  (m >= n)
    !
    ! OUT :  rinvqt = R^{-1}*Q^t, which gives the solution as x = R^{-1}*Q^t*b.
    !
    !*****************************************************************************
    subroutine qr_factorization(a,rinvqt,m,n)

        implicit none
    
        !Input
        integer ,                 intent( in) :: m, n
        real(p2), dimension(m,n), intent( in) :: a
    
        !Output
        real(p2), dimension(n,m), intent(out) :: rinvqt
            
        ! Local variables
        !
        ! Note: Think if you can reduce the number of
        !       variables below to save memory.
        
        integer                  :: i, j, k, ii
    
        real(p2), dimension(:,:), pointer :: r
        real(p2)                          :: abs_rk, sign_rkk, wTw
        real(p2), dimension(m)            :: w, rk
        real(p2), dimension(m,m)          :: qt, wwT
    
        real(p2), dimension(:,:), pointer ::  r_nxn
        real(p2), dimension(:),   pointer   ::   y, b
        real(p2)                 ::    rhs
        integer,  dimension(n)   :: zeroCheck
        integer                  :: nonZeroColumns
        real(p2), dimension(n,m) :: rinvqt_intermediate
        
        do i = 1,n
            zeroCheck(i) = 0
        end do
        nonZeroColumns = 0
        do i = 1,n
            if (.not.(vector_norm(a(:,i),m)) == 0) then
                zeroCheck(i) = 1
                nonZeroColumns = nonZeroColumns + 1
            end if
        end do
        if (nonZeroColumns == 0) then
            write(*,*) " error dx = dy = dz = 0" ! This really shouldn't happen...
        end if
        allocate(r(m,nonZeroColumns))
        ii = 1
        do i = 1,n
            if (zeroCheck(i) == 1) then
                r(:,ii) = a(:,i)
                ii = ii + 1
            end if
        end do

        if (m < n) then
        write(*,*) " Underdetermined system detected... m < n: "
        write(*,*) "   m =  ", m
        write(*,*) "   n =  ", n
        write(*,*) " qr_factorization() not designed to solve such a problem... Stop. "
        stop
        endif
    
        !-------------------------------------------------------
        ! Initialization: R = A
    
        ! r = a ! not anymore
    
        !-------------------------------------------------------
        ! Initialization: Qt = I
        
        qt = zero
    
        do i = 1, m
            qt(i,i) = one
        end do
    
        !-------------------------------------------------------
        ! Apply reflection to each column of R, and generate
        ! the final upper triangular matrix R and the transpose
        ! Qt of the orthonormal matrix Q.
    
        column_loop : do k = 1, nonZeroColumns
    
            !Our target are the elements below the (k,k) element
            !in k-th column, i.e., r(k:m).
            !So, rk gets shorter as we move on (as k increases).
    
            rk      = zero
            rk(k:m) = r(k:m,k)
    
            !Reflector Hk will zero out all the elements below r(k).
        
            !Compute the length of rk and the sign of the kth element.
    
              abs_rk = sqrt( dot_product(rk,rk) )
            sign_rkk = sign( one, rk(k) )
    
            !Define the reflecting vector w:   w = |rk|*(1,0,0,...,0)-rk
            !                               or w =-|rk|*(1,0,0,...,0)-rk
            !We switch the reflection (there are two possible ones)
            !to avoid w = 0 that can happen if rk=(1,0,0,...,0).
    
            w      = zero
            w(k)   = -sign_rkk*abs_rk
            w(k:m) = w(k:m) - rk(k:m)
    
            !Compute the length^2 of w: wt*w = [x,x,...,x]|x| = dot product = scalar.
            !                                             |x|
            !                                             |.|
            !                                             |.|
            !                                             |x|
        
            wTw = dot_product(w,w)
        
            !Compute the dyad of w: w*wt = |x|[x,x,...,x] = mxm matrix.
            !                              |x|
            !                              |.|
            !                              |.|
            !                              |x|
        
            do i = 1, m
                do j = 1, m
                    wwT(i,j) = w(i)*w(j)
                end do
            end do
        
            !We now apply the reflector matrix Hk = I-2*wwt/wTw,
            !and update R and Qt.
        
            !Update  R:  R = Hk*R  = (I-2*wwt/wTw)*R  = R-2*(wwt*R)/wTw
        
            r  =  r - two*matmul(wwT,r)/wTw
        
            !Update Qt: Qt = Hk*Qt = (I-2*wwt/wTw)*Qt = Qt-2*(wwt*Qt)/wTw
        
            qt = qt - two*matmul(wwT,qt)/wTw
        
        end do column_loop
    
        !-------------------------------------------------------
        ! Compute rinvqt(1:n,1:m) = R_{nxn}^{-1} * Q_{nxm}^t by
        ! solving R_{nxn} * rinvqt(1:n,k) = Q_{nxm}^t(1:n,k)
        ! for k=1,n. We can solve it easily by back substitution
        ! since R_{nxn} is upper triangular.
        
        allocate(r_nxn(nonZeroColumns,nonZeroColumns))
        r_nxn =  r(1:nonZeroColumns,1:nonZeroColumns)
    
        allocate(y(nonZeroColumns))
        allocate(b(nonZeroColumns))
        rinvqt = zero ! initialize as zero since some values won't get assigned later
        rinvqt_intermediate = zero
        do k = 1, m
    
            !Solve r*y = b, where y is the k-th column of rinvqt.
    
            b = qt(1:nonZeroColumns,k)
    
            !Solve the lower right equation.
    
            rhs = b(nonZeroColumns)
            y(nonZeroColumns) = rhs/r_nxn(nonZeroColumns,nonZeroColumns)
    
            !Go up and solve.
            do i = nonZeroColumns-1, 1, -1
    
                !Take all known parts (j=i+1,n) to the rhs.
        
                !RHS is known, of course.
                rhs = b(i)
                !Below are all known since the solutions y(j=i+1,n) has already been solved.
                do j = i+1, nonZeroColumns
                    rhs = rhs - r_nxn(i,j)*y(j)
                end do
        
                !Divide the rhs by the coefficient of the (i,i) part.
                y(i) = rhs/r_nxn(i,i)
        
            end do
    
            !The soluton x is the k-th column of rinvqt.
            rinvqt_intermediate(:nonZeroColumns,k) = y(:)    
        end do
        nonZeroColumns = 0
        do i = 1,n
            if (zeroCheck(i) == 1) then
                nonZeroColumns = nonZeroColumns + 1
                rinvqt(i,:) = rinvqt_intermediate(nonZeroColumns,:)
            end if
        end do   
        deallocate(r,r_nxn,y,b) 
    end subroutine qr_factorization

        ! Function for computing the L2 vector norm from http://me.rice.edu/~akin/OOP_Copyrighted/4_Features_of_Lang/vector_norm.f90
    ! modified from a standalone program
    function vector_norm(x, n)
        !   A simple vector norm program, Fortran 90 version
        use common   , only : p2
        implicit none
            !integer, parameter :: dp = selected_real_kind(14) ! 14 digits
            real(p2), dimension(:), intent(in)  :: x    
            integer, intent(in)                 :: n

            real(p2) :: vector_norm

            vector_norm = sqrt ( sum ( x(:n)*x(:n) ))   ! L2 norm

           
    end function vector_norm

end module direct_solve