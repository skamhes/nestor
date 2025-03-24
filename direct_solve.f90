module direct_solve

    use common , only : p2, zero, one, two

    implicit none

    public :: qr_factorization      ! QR factorization using Householder reflections 
    ! https://en.wikipedia.org/wiki/QR_decomposition#Using_Householder_reflections:~:text=Using-,Householder,-reflections%5Bedit 
    public :: gewp_solve            ! Gauss elimination for inverting diagonal blocks
    public :: safe_invert_scalar

    private

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


    !****************************************************************************
    !* ------------------ GAUSS ELIMINATION WITH PIVOTING ---------------------
    !*
    !*  This computes the inverse of an (nm)x(nm) matrix "ai" and also
    !*  computes the solution to a given lienar system.
    !*
    !*  IN :       ai = An (nm)x(nm) matrix whoise inverse is sought.
    !*             bi = A vector of (nm): Right hand side of the linear sytem
    !*             nm = The size of the matrix "ai"
    !*
    !* OUT :
    !*            sol = Solution to the linear system: ai*sol=bi
    !*        inverse = the inverse of "ai".
    !*       idetstat = 0 -> inverse successfully computed
    !*                  1 -> THE INVERSE DOES NOT EXIST (det=0).
    !*                  2 -> No unique solutions exist.
    !*****************************************************************************
    subroutine gewp_solve(ai,nm, inverse,idetstat)

        implicit none
      
        integer , parameter ::    p2 = selected_real_kind(15) ! Double precision
        real(p2), parameter ::  zero = 0.0_p2
        real(p2), parameter ::   one = 1.0_p2
      
        integer ,                   intent( in) :: nm
        real(p2), dimension(nm,nm), intent( in) :: ai
      
        real(p2), dimension(nm,nm), intent(out) :: inverse
        integer ,                   intent(out) :: idetstat
      
        real(p2), dimension(nm,nm+1) :: a
        real(p2), dimension(nm)      :: x
        integer , dimension(nm)      :: nrow
        integer                      :: I,J,K,pp,m
      
        do m = 1, nm
            !*****************************************************************************
            !* Set up the matrix a
            !*****************************************************************************
            
            do J=1,nm
                do I=1,nm
                a(I,J) = ai(I,J)
                end do
            end do
        
            do k=1,nm
                a(k,nm+1)=zero; nrow(k)=k
            end do
            a(m,nm+1)=one
        
            !*****************************************************************************
            !* HONA IKOKA..... 
            !*****************************************************************************
            do j=1,nm-1
            !*****************************************************************************
            !* FIND SMALLEST pp FOR a(pp,j) IS MAXIMUM IN JTH COLUMN.
            !***************************************************************************** 
                call findmax(nm,j,pp,a,nrow)
                !*****************************************************************************
                !* IF a(nrow(p),j) IS zero, THERE'S NO UNIQUE SOLUTIONS      
                !*****************************************************************************
                if (abs(a(nrow(pp),j)) < epsilon(one)) then
                    write(6,*) 'THE INVERSE DOES NOT EXIST.'
                    idetstat = 1
                    return
                endif
                !*****************************************************************************
                !* IF THE MAX IS NOT A DIAGONAL ELEMENT, SWITCH THOSE ROWS       
                !*****************************************************************************
                if (nrow(pp) .ne. nrow(j)) then
                    call switch(nm,j,pp,nrow)
                else
                endif  
                !*****************************************************************************
                !* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
                !***************************************************************************** 
                call eliminate_below(nm,j,a,nrow)
        
            end do
            !*****************************************************************************
            !* CHECK IF a(nrow(N),N)=0.0 .
            !*****************************************************************************
            if (abs(a(nrow(nm),nm)) < epsilon(one)) then
                write(6,*) 'NO UNIQUE SOLUTION EXISTS!'
                idetstat = 2
                return
            else
            endif
            !*****************************************************************************
            !* BACKSUBSTITUTION!
            !*****************************************************************************
            call backsub(nm,x,a,nrow)
            !*****************************************************************************
            !* STORE THE SOLUTIONS, YOU KNOW THEY ARE INVERSE(i,m) i=1...
            !*****************************************************************************
            do i=1,nm
                inverse(i,m)=x(i)
            end do
            !*****************************************************************************
        end do
      
        idetstat = 0
      
        return
      
        !*****************************************************************************
    end subroutine gewp_solve
      
      !*****************************************************************************
      !* Four subroutines below are used in gewp_solve() above.
      !*****************************************************************************
      !* FIND MAXIMUM ELEMENT IN jth COLUMN 
      !***************************************************************************** 
            subroutine findmax(nm,j,pp,a,nrow)
      
            implicit none
      
            integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
            integer , intent( in) :: nm
            real(p2), intent( in) :: a(nm,nm+1)
            integer , intent( in) :: j,nrow(nm)
            integer , intent(out) :: pp
            real(p2)              :: max
            integer               :: i
      
                  max=abs(a(nrow(j),j)); pp=j
      
                 do i=j+1,nm
      
                   if (max < abs(a(nrow(i),j))) then
      
                        pp=i; max=abs(a(nrow(i),j))
      
                   endif
      
                 end do
      
            return
      
            end subroutine findmax
      !*****************************************************************************
      !* SWITCH THOSE ROWS       
      !*****************************************************************************
            subroutine switch(nm,j,pp,nrow)
      
            implicit none
      
            integer, intent(   in) :: nm,j,pp
            integer, intent(inout) :: nrow(nm)
            integer                :: ncopy
      
            if (nrow(pp).ne.nrow(j)) then
      
               ncopy=nrow(j)
               nrow(j)=nrow(pp)
               nrow(pp)=ncopy
      
            endif
      
            return
      
            end subroutine switch
      !*****************************************************************************
      !* ELIMINATE ALL THE ENTRIES BELOW THE DIAGONAL ONE
      !*(Give me j, the column you are working on now)
      !***************************************************************************** 
            subroutine eliminate_below(nm,j,a,nrow)
      
            implicit none
      
            integer , parameter     :: p2 = selected_real_kind(15) ! Double precision
            real(p2), parameter     :: zero = 0.0_p2
            integer , intent(   in) :: nm
            real(p2), intent(inout) :: a(nm,nm+1)
            integer , intent(   in) :: j,nrow(nm)
            real(p2)                :: m
            integer                 :: k,i
      
            do i=j+1,nm
      
              m=a(nrow(i),j)/a(nrow(j),j)
              a(nrow(i),j)=zero
      
                do k=j+1,nm+1
                  a(nrow(i),k)=a(nrow(i),k)-m*a(nrow(j),k)
                end do
      
            end do
      
            return
      
            end subroutine eliminate_below
      !*****************************************************************************
      !* BACKSUBSTITUTION!
      !*****************************************************************************
            subroutine backsub(nm,x,a,nrow)
      
            implicit none
      
            integer , parameter   :: p2 = selected_real_kind(15) ! Double precision
            real(p2), parameter   :: zero = 0.0_p2
      
            integer , intent( in) :: nm
            real(p2), intent( in) :: a(nm,nm+1)
            integer , intent( in) :: nrow(nm)
            real(p2), intent(out) :: x(nm)
            real(p2)              :: sum
            integer               :: i,k
      
            x(nm)=a(nrow(nm),nm+1)/a(nrow(nm),nm)
      
            do i=nm-1,1,-1
      
               sum=zero
      
                 do k=i+1,nm
      
                    sum=sum+a(nrow(i),k)*x(k)
      
                 end do
      
            x(i)=(a(nrow(i),nm+1)-sum)/a(nrow(i),i)
      
            end do
      
            return
      
            end subroutine backsub
      !*********************************************************************

    pure elemental function safe_invert_scalar(scal) result(inv)
        use common , only : p2, one

        ! This function prevents divide by zero errors by capping very small inputs at a set cutoff from zero (either positive or
        ! negative).
        ! note this function assumes the input is not NaN or -NaN
        
        implicit none
        real(p2), intent(in)    :: scal
        real(p2)                :: inv

        real(p2), parameter :: cutoff = 1.0e-012
        
        ! This method avoids any branching.  I don't know if it's actually faster but it looks faster and that's what really counts
        inv = one / ( sign(one,scal) * max( abs(scal) , cutoff ) )
    end function safe_invert_scalar
end module direct_solve