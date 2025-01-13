module res_sa

    use common , only : p2

    private

    public sa_invFlux
    public sa_viscFlux
    public sa_sourceTerms

    contains

    subroutine sa_invFlux(mu1, mu2, q1, q2, gradmu1, gradmu2, n12, xc1, yc1, zc1, xc2, yc2, zc2, &
                             xm, ym, zm, phi1, phi2, mu_flux, jac1, jac2 )

        use config , only : rans_accuracy

        use common , only : half, zero

        use solution , only : nq
                            
        implicit none

        real(p2),               intent(in) :: mu1, mu2
        real(p2), dimension(:), intent(in) :: q1,  q2
        real(p2), dimension(:), intent(in) :: gradmu1, gradmu2
        real(p2), dimension(:), intent(in) :: n12
        real(p2),               intent(in) :: xc1, yc1, zc1, xc2, yc2, zc2
        real(p2),               intent(in) :: xm, ym, zm
        real(p2),               intent(in) :: phi1, phi2

        real(p2),               intent(out):: mu_flux
        real(p2),               intent(out):: jac1, jac2

        real(p2) :: muL, muR

        real(p2), dimension(nq) :: qF
        real(p2)                :: vL, vR, vF


        integer, parameter :: iu = 2
        integer, parameter :: iv = 3
        integer, parameter :: iw = 4

        if (rans_accuracy == 2) then
            muL = mu1 + phi1 * ( gradmu1(1)*(xm-xc1) + gradmu1(2)*(ym-yc1) + gradmu1(3)*(zm-zc1) ) ! gradmu <=> gradmu (var) 
            muR = mu2 + phi2 * ( gradmu2(1)*(xm-xc2) + gradmu2(2)*(ym-yc2) + gradmu2(3)*(zm-zc2) ) ! u <=> mu (vars)
        else
            muL = mu1 
            muR = mu2
        end if

        qF(iu) = half * (q1(iu) + q2(iu))
        qF(iv) = half * (q1(iv) + q2(iv))
        qF(iw) = half * (q1(iw) + q2(iw))

        vF = qF(iu) * n12(1) + qF(iv) * n12(2) + qF(iw) * n12(3)
        vL = max(vF,zero)
        vR = min(vF,zero)

        mu_flux = vL * muL + vR * muR
        jac1 = vL
        jac2 = vR

    end subroutine sa_invFlux
end module res_sa