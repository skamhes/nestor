module viscosity

    use common , only : p2
    
    private

    public compute_viscosity, compute_viscosity_ddt

    contains 



    pure elemental function compute_viscosity(T) result(mu)

        use common , only : three_half, one

        use config , only : M_inf, Re_inf, sutherland_constant, reference_temp
        
        use solution_vars , only : C0, mre

        implicit none

        real(p2), intent(in) :: T
        real(p2)             :: mu

        ! note: C0 = C0 / T_inf in reality.  SEE: set_initial_solution in initialize.f90
        mu =  mre * (one + C0) / (T + C0)*T**(three_half)

    end function compute_viscosity

    pure elemental function compute_viscosity_ddt(T) result(mu)

        use common , only : three_half, one

        use config , only : M_inf, Re_inf, sutherland_constant, reference_temp
        
        use solution_vars , only : C0, mre

        use ad_operators

        implicit none

        type(derivative_data_type_df5), intent(in) :: T
        type(derivative_data_type_df5)             :: mu

        ! note: C0 = C0 / T_inf in reality
        mu =  mre * (one + C0) / (T + C0)*T**(three_half)

    end function compute_viscosity_ddt
end module viscosity