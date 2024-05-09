module recfast_wrapper
    use iso_c_binding, only: c_double, c_int
    use recfast_module, only: recfast_func
    implicit none

    contains

        subroutine recfast_c(Omega_b, Omega_c, Omega_L, H0, T_CMB, Yp, H_switch, He_switch, &
                             z_initial, z_final, tol, Nz, z_array, x_array) bind(c)
        real(c_double), intent(in) :: Omega_b
        real(c_double), intent(in) :: Omega_c
        real(c_double), intent(in) :: Omega_L
        real(c_double), intent(in) :: H0
        real(c_double), intent(in) :: T_CMB
        real(c_double), intent(in) :: Yp
        integer(c_int), intent(in) :: H_switch
        integer(c_int), intent(in) :: He_switch
        real(c_double), intent(in) :: z_initial
        real(c_double), intent(in) :: z_final
        real(c_double), intent(in) :: tol
        integer(c_int), intent(in) :: Nz
        real(c_double), intent(out) :: z_array(Nz)
        real(c_double), intent(out) :: x_array(Nz)
        call recfast_func(Omega_b, Omega_c, Omega_L, H0, T_CMB, Yp, H_switch, He_switch, &
                          z_initial, z_final, tol, Nz, z_array, x_array)
    end subroutine recfast_c
end module recfast_wrapper

