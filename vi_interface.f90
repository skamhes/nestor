module vi_interface

    use iso_c_binding, only: c_double, c_int

    implicit none

    public intrinsic_grad
    
    interface
        subroutine old_intrinsic_grad(q_, ci_, ck_, cxyz_, ccgradq_) bind(C, name="old_intrinsic_grad")
            use iso_c_binding, only: c_double, c_int
            real(c_double), dimension(*) :: q_ 
            real(c_double), dimension(*) :: cxyz_
            real(c_double), dimension(*) :: ccgradq_ ! Passed by reference
            integer(c_int), value          :: ci_, ck_

        end subroutine old_intrinsic_grad
    end interface


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

    interface
        subroutine transpose_3x5(cT__,c__) bind(C, name="transpose_3x5")
            use iso_c_binding, only: c_double, c_int
            real(c_double), dimension(*) :: cT__ 
            real(c_double), dimension(*) :: c__
        end subroutine transpose_3x5
    end interface
end module vi_interface