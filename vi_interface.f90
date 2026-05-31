module vi_interface

    use iso_c_binding, only: c_double, c_int

    implicit none

    public intrinsic_grad
    
    interface
        subroutine intrinsic_grad(q_, ci_, ck_, cxyz_, ccgradq_) bind(C, name="intrinsic_grad")
            use iso_c_binding, only: c_double, c_int
            real(c_double), dimension(*) :: q_ 
            real(c_double), dimension(*) :: cxyz_
            real(c_double), dimension(*) :: ccgradq_ ! Passed by reference
            integer(c_int), value          :: ci_, ck_

        end subroutine intrinsic_grad
    end interface

end module vi_interface