module vi_interface

    use iso_c_binding, only: c_double, c_int

    implicit none

    public intrinsic_grad

    interface
        subroutine intrinsic_grad(q_, icell_, nnghbrs_, ckn_, cf_, & ! Internal cells
                                  qg_,        nbf_,           gcf_, & ! boundary cells
                                  ccgradq_) bind(C, name="intrinsic_grad") ! OUTPUT
            use iso_c_binding, only: c_double, c_int, c_ptr
            real(c_double), dimension(*) :: q_
            type(c_ptr),    dimension(*) :: qg_
            integer(c_int), value        :: icell_, nnghbrs_, nbf_
            integer(c_int), dimension(*) :: ckn_
            real(c_double), dimension(*) :: cf_, gcf_
            real(c_double), dimension(*) :: ccgradq_ ! Passed by reference
        end subroutine intrinsic_grad
    end interface

end module vi_interface