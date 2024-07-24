module gauss_seidel

    implicit none

    public :: FORWARD, BACKWARD
    public :: gauss_seidel_sweep
    integer, parameter :: FORWARD  = 1
    integer, parameter :: BACKWARD = 1

    contains

    subroutine gauss_seidel_sweep(nq,res,V,R,C,Dinv,nnz,omega_lrelax,correction, linear_res_norm)

        use common          , only : p2, zero

        use grid            , only : ncells
        
        implicit none

        !INPUT
        integer         :: nq ! Number of equations
        integer ,                           intent(in)   :: nnz
        real(p2), dimension(   nq,nnz),     intent(in)   :: res    ! RHS (b)
        real(p2), dimension(nq,nq,nnz),     intent(in)   :: V    ! Values of A
        integer , dimension(nnz      ),     intent(in)   :: C    ! Column index of A
        integer , dimension(ncells+1 ),     intent(in)   :: R    ! Start index of A
        real(p2), dimension(nq,nq,ncells),  intent(in)   :: Dinv ! Inverse of A(i,i)
        real(p2),                           intent(in)   :: omega_lrelax
        !INOUT
        real(p2), dimension(nq,ncells ), intent(inout)   :: correction
        !OUTPUT
        real(p2), dimension(nq)        , intent(out)     :: linear_res_norm
        
        ! local
        real(p2), dimension(nq) :: b ! rhs
        integer :: i,k
        real(p2), dimension(nq) :: linear_res

        linear_res_norm = zero

        gs_loop : do i = 1,ncells ! loop through rows
            ! Form the right hand side of GS: [ sum( off_diagonal_block*du ) - residual ]
            b = -res(:,i)
            gs_row_loop : do k = R(i),(R(i+1)-1)
                ! Add RHS from off diagonal terms and du (du = zero to start and will be updated as we go)
                if ( C(k) .NE. i) then
                    b = b - matmul(V(:,:,k),correction(:,C(k)))
                end if
            end do gs_row_loop
            ! ! Update du by the GS relaxation:
                    !
                    ! e.g., for 3 nghbrs, perform the relaxation in the form:
                    !
                    !                     diagonal block        sum of off-diagonal block contributions
                    !       dUj = omega*{ [V/dtj+dR/dUj]^{-1}*(-[dRj/dU1]*dU1 -[dRj/dU2]*dU2 -[dRj/dU3]*dU3 -Res_j) - dUj }
            linear_res = matmul(Dinv(:,:,i), b) - correction(:,i)
            correction(:,i) = correction(:,i) + omega_lrelax * linear_res
            linear_res_norm(:) = linear_res_norm(:) + abs(linear_res)
        end do gs_loop

        !---------------------------------------------------------
        linear_res_norm(:) = linear_res_norm(:) / real(ncells, p2)

    end subroutine gauss_seidel_sweep
end module gauss_seidel