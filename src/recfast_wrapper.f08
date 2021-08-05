module recfast_wrapper
    use iso_c_binding, only: c_double, c_int
    use recfast_module, only: recfast_func
    implicit none

    contains

        subroutine c_recfast(OmegaB, OmegaC, OmegaL, H0inp, Tnow, Yp, Hswitch_in, Heswitch_in, &
                             zinitial, zfinal, tol, Nz, z_array, x_array) bind(c)
        integer(c_int), intent(in) :: Nz
        integer(c_int), intent(in) :: Hswitch_in
        integer(c_int), intent(in) :: Heswitch_in
        real(c_double), intent(in) :: OmegaB
        real(c_double), intent(in) :: OmegaC
        real(c_double), intent(in) :: OmegaL
        real(c_double), intent(in) :: H0inp
        real(c_double), intent(in) :: Tnow
        real(c_double), intent(in) :: Yp
        real(c_double), intent(in) :: zinitial
        real(c_double), intent(in) :: zfinal
        real(c_double), intent(in) :: tol
        real(c_double), intent(out) :: z_array(Nz)
        real(c_double), intent(out) :: x_array(Nz)
        call recfast_func(OmegaB, OmegaC, OmegaL, H0inp, Tnow, Yp, Hswitch_in, Heswitch_in, &
                zinitial, zfinal, tol, Nz, z_array, x_array)
    end subroutine c_recfast
end module recfast_wrapper

