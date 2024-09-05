module gcr

    ! Method for computing the generalized conjugate residual, a Jacobian Free Newton-Krylov method.
    ! This is based on the paper: https://doi.org/10.2514/6.2021-0857 (A free version can be found on the NASA TRS).

    implicit none

    private 

    public :: run_gcr

    contains

    subroutine run_gcr()

        use common      , only : p2, zero, one, half

        use config      , only : gcr_max_projections

        use grid        , only : ncells, cell

        use solution    , only : res, nq, ndim, jac, Q, inv_ncells, gamma, gammamo, gmoinv, dtau

        use jacobian    , only : compute_jacobian

        use linear_solver, only: linear_relaxation

        real(p2), dimension(nq,ncells,gcr_max_projections) :: gcr_precond_correction
        real(p2), dimension(nq,ncells)                     :: gcr_final_correction
        real(p2), dimension(nq,ncells,gcr_max_projections) :: gcr_search_direction

        real(p2), dimension(gcr_max_projections)           :: gamma_projection ! gamma_k from EQ. 14
        real(p2), dimension(gcr_max_projections)           :: gcr_search_direction_mag
        real(p2)                                           :: gcr_inner_prod
        real(p2), dimension(nq,ncells)                     :: Q0    ! Initial primative vector at the start of outer iteration
        real(p2), dimension(nq,ncells)                     :: Rprev ! Residual from the previous search direction
        real(p2)                                           :: correction_rms
        real(p2)                                           :: correction_norm ! l2 norm of preconditioned correction
        real(p2)                                           :: eps_frechet
        
        real(p2), dimension(5,5)    :: preconditioner
        real(p2)                    :: theta, rho_p, rho_T, rho
        real(p2)                    :: H, alpha, beta, lambda, absu, UR2inv


        integer :: kdir, jdir ! projection direction indices
        integer :: icell
        integer :: ivar
        integer :: ii

        ! Store the initial values.  For now this isn't terribly efficient.  Maybe I'll change it later...
        Rprev = - res
        Q0    = Q


        ! Generate the initial Approximate Jacobian and compute an initial update (preconditioner)
        call compute_jacobian
        
        project_loop : do kdir = 1,gcr_max_projections
            call linear_relaxation(nq, jac, res, gcr_precond_correction(:,:,kdir))

            ! Compute Frechet Derivative
            correction_rms = zero
            frechet : do icell = 1,ncells
                do ivar = 1,nq
                    correction_rms = correction_rms + gcr_precond_correction(ivar,icell,kdir)**2
                    correction_norm = correction_norm + gcr_precond_correction(ivar,icell,kdir)
                end do
            end do frechet
            correction_rms = sqrt( correction_rms*inv_ncells )
            correction_norm = correction_norm * inv_ncells
            eps_frechet = max(one,correction_rms)

            q(:,:) = q(:,:) + eps_frechet * gcr_precond_correction(:,:,kdir) / correction_norm

            call compute_residual

            ! Generate new search direction

            do icell = 1,ncells
                H = ((q(5,icell))**2)*gmoinv + half * ( q(2,icell)**2 + q(3,icell)**2 + q(4,icell)**2 )
                rho_p = gamma/q(5,icell)
                rho_T = - (q(1,icell)*gamma)/(q(5,icell)**2)
                rho = q(1,icell)*gamma/q(5,icell)
                UR2inv = one ! will be 1/uR2(i)
                theta = (UR2inv) - rho_T*(gammamo)/(rho)
                
                ! Note transposing this assignment would likely be marginally faster if slightly less easy to read
                preconditioner(1,:) = (/ theta,            zero,           zero,           zero,           rho_T                  /)
                preconditioner(2,:) = (/ theta*q(2,icell), rho,            zero,           zero,           rho_T*q(2,icell)       /)
                preconditioner(3,:) = (/ theta*q(3,icell), zero,           rho,            zero,           rho_T*q(3,icell)       /)
                preconditioner(4,:) = (/ theta*q(4,icell), zero,           zero,           rho,            rho_T*q(4,icell)       /)
                preconditioner(5,:) = (/ theta*H-one,      rho*q(2,icell), rho*q(3,icell), rho*q(4,icell), rho_T*H+rho/(gamma-one)/)
    
                do ii = 1,nq
                    preconditioner(ii,ii) = preconditioner(ii,ii) * cell(icell)%vol/dtau(icell)
                end do
            end do

        enddo project_loop


    end subroutine

end module gcr