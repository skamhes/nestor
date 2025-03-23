module viscosity

    use common , only : p2
    
    private

    public compute_viscosity, compute_viscosity_ddt
    public C0, mu_norm
    
    real(p2) :: C0
    real(p2) :: mu_norm

    contains 



    pure elemental function compute_viscosity(T) result(mu)

        use common , only : three_half, one

        use config , only : M_inf, Re_inf, sutherland_constant, reference_temp

        implicit none

        real(p2), intent(in) :: T
        real(p2)             :: mu

        ! note: C0 = C0 / T_inf in reality
        mu =  mu_norm * (one + C0) / (T + C0)*T**(three_half)

    end function compute_viscosity

    pure elemental function compute_viscosity_ddt(T) result(mu)

        use common , only : three_half, one

        use config , only : M_inf, Re_inf, sutherland_constant, reference_temp

        use ad_operators

        implicit none

        type(derivative_data_type_df5), intent(in) :: T
        type(derivative_data_type_df5)             :: mu

        ! note: C0 = C0 / T_inf in reality
        mu =  mu_norm * (one + C0) / (T + C0)*T**(three_half)

    end function compute_viscosity_ddt
end module viscosity