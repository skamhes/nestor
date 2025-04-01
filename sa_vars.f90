module sa_vars

    use common , only : p2, one, two

    implicit none

    real(p2), parameter :: cb1    = 0.1355_p2
    real(p2), parameter :: cb2    = 0.622_p2
    real(p2), parameter :: KAPPA  = 0.41_p2
    real(p2), parameter :: cw2    = 0.3_p2
    real(p2), parameter :: cw3    = 2.0_p2
    real(p2), parameter :: cv1    = 7.1_p2
    real(p2), parameter :: ct3    = 1.2_p2
    real(p2), parameter :: ct4    = 0.5_p2
    real(p2), parameter :: SIGMA  = 2.0_p2 / 3.0_p2 ! 1/sigma, sigma = 2/3
    real(p2), parameter :: iSIGMA = 1.0_p2 / SIGMA ! 1/sigma, sigma = 2/3
    real(p2), parameter :: c2     = 0.7_p2
    real(p2), parameter :: c22    = 0.49_p2
    real(p2), parameter :: c3     = 0.9_p2

    ! Computed
    real(p2), parameter :: cw1    = ( cb1 / KAPPA**2 ) + ( ( one + cb2 ) * iSIGMA )
    real(p2), parameter :: c3m2c2 = c3 - two * c2
    real(p2), parameter :: cw36   = cw3**6
    real(p2), parameter :: cv13   = cv1**3

end module sa_vars