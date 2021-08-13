!##############################################################################
! Integrator for Cosmic Recombination of Hydrogen and Helium,
! developed by Douglas Scott (dscott@astro.ubc.ca)
! based on calculations in the papers Seager, Sasselov & Scott
! (ApJ, 523, L1, 1999; ApJS, 128, 407, 2000)
! and "fudge" updates in Wong, Moss & Scott (2008).
!
! Permission to use, copy, modify and distribute without fee or royalty at
! any tier, this software and its documentation, for any purpose and without
! fee or royalty is hereby granted, provided that you agree to comply with
! the following copyright notice and statements, including the disclaimer,
! and that the same appear on ALL copies of the software and documentation,
! including modifications that you make for internal use or for distribution:
!
! Copyright 1999-2010 by University of British Columbia.  All rights reserved.
!
! THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO
! REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
! BY WAY OF EXAMPLE, BUT NOT LIMITATION,
! U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF
! MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
! THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
! ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
!
!##############################################################################
!
!N  Name: RECFAST
!V  Version: 1.8.0
!
!P  Purpose:
!P      Calculate ionised fraction as a function of redshift.
!P      Solves for H and He simultaneously, and includes
!P      H "fudge factor" for low z effect, as well as
!P      HeI fudge factor.
!
!D  Description:
!D      Solves for ionisation history since recombination
!D      using the equations in Seager, Sasselov & Scott (ApJ, 1999).
!D      The Cosmological model can be flat or open.
!D      The matter temperature is also followed, with an update from
!D      Moss & Scott (2009).
!D      The values for \alpha_B for H are from Hummer (1994).
!D      The singlet HeI coefficient is a fit from the full code.
!D      Additional He "fudge factors" are as described in Wong, Moss
!D      and Scott (2008).
!D      Extra fitting function included (in optical depth) to account
!D      for extra H physics described in Rubino-Martin et al. (2010).
!D      Care is taken to use the most accurate constants.
!D      Note that some parameters are fixed (e.g. N_nu=3, nu's are
!D      massless, w=-1, etc.) - some users may want to explictly
!D      input their own H(z) to account for extra physics.
!D      This is provided as a PROGRAM, which can be easily converted
!D      to a SUBROUTINE for use in CMB Boltzmann codes.
!
!F  File & device access:
!F      Unit    /I, IO, O /Name (if known)
!
!M  Modules used:
!M      precision (module to declare `dp` for double precision)
!M      ode_solver (module containing the ODE solver `dverk`)
!M      constants (module for mathematical and physical constants and chemical data)
!M      fudgefit (module for fudge factors and fitting parameters)
!M      cosmo_input  (module to set up cosmological input parameters and make them available throughout)
!M      recombination (main recfast module)
!M
!M  Subroutines called:
!M      dverk (numerical integrator)
!M      get_init (initial values for ionization fractions)
!M      ion (ionization and Temp derivatives)
!
!H  History:
!H  CREATED     (simplest version) 19th March 1989
!H  RECREATED   11th January 1995
!H          includes variable Cosmology
!H          uses dverk integrator
!H          initial conditions are Saha
!H  TESTED      a bunch, well, OK, not really
!H  MODIFIED
!H      Jan 1995 (include Hummer's 1994 alpha table)
!H      Jan 1995 (include new value for 2s-1s rate)
!H      Jan 1995 (expand comments)
!H      Mar 1995 (add Saha for Helium)
!H      Aug 1997 (add HeII alpha table)
!H      Jul 1998 (include OmegaM correction and H fudge factor)
!H      Nov 1998 (change Trad to Tmat in Rup)
!H      Jan 1999 (tidied up for public consumption)
!H      Sep 1999 (switch to formula for alpha's, fix glitch)
!H      Feb 2000 (fixed overflow problem in He_Boltz)
!H      Oct 2001 (fixed OmegaM in z_eq)
!H      Jun 2003 (fixed error in Rdown_He formula)
!H      Jun 2003 (fixed He recombination coefficients)
!H      Jun 2003 (comments to point out fixed N_nu etc.)
!H      Oct 2006 (included new value for G)
!H      Oct 2006 (improved m_He/m_1H to be "not4")
!H      Oct 2006 (fixed error, x for x_H in part of f(1))
!H      Jan 2008 (improved HeI recombination effects, including HeI rec. fudge factor)
!H      Feb 2008 (avoid calculating gamma_2Ps and gamma_2Pt when x_H close to 1.0)
!H      Aug 2008 (correction for x_H when Heflag=2 and Helfag>=5 to be smoother)
!H      Sep 2008 (added extra term to make transition smoother for Tmat evolution)
!H      Jan 2010 (added fitting function to modify K to match x_e(z) for new H physics)
!H      Jul 2012 (modified fudge factors for better match to codes with more detailed physics)
!H      Sep 2012 (fixed "fu" at low z to match modifications)
!H      Jul 2021 (update constants, convert to Modern Fortran (2008) and make Python wrapper)
!
!##############################################################################


!##############################################################################
module fudgefit
    use precision, only : dp
    implicit none

    ! Gaussian fits for extra H physics (fit by Adam Moss, modified by Antony Lewis):
    real(dp), parameter :: AGauss1 = -0.14_dp  ! Amplitude of 1st Gaussian for the H fudging
    real(dp), parameter :: AGauss2 = 0.079_dp  ! Amplitude of 2nd Gaussian for the H fudging
    real(dp), parameter :: zGauss1 = 7.28_dp   ! ln(1 + z) central value of 1st Gaussian
    real(dp), parameter :: zGauss2 = 6.73_dp   ! ln(1 + z) central value of 2nd Gaussian
    real(dp), parameter :: wGauss1 = 0.18_dp   ! Width of 1st Gaussian
    real(dp), parameter :: wGauss2 = 0.33_dp   ! Width of 2nd Gaussian

    ! Matter departs from radiation when t(Th) > H_frac * t(H)
    ! choose some safely small number:
    real(dp), parameter :: H_frac = 1.e-3_dp   ! H_frac: follow Tmat when t_Compton / t_Hubble > H_frac

    integer  :: Hswitch   ! Hswitch is an integer for modifying the H recombination
    integer  :: Heswitch  ! Heswitch is an integer for modifying HeI recombination
    real(dp) :: fu        ! fu is a "fudge factor" for H, to approximate low z behaviour
    real(dp) :: b_He      ! b_He is a "fudge factor" for HeI, to approximate higher z behaviour

    contains

    subroutine set_switch(Hswitch_in, Heswitch_in)
        integer, intent(in) :: Hswitch_in
        integer, intent(in) :: Heswitch_in

        Hswitch = Hswitch_in
        Heswitch = Heswitch_in

        ! Fudge factor to approximate the low z out of equilibrium effect
        if (Hswitch == 0) then
            fu = 1.14_dp
        else
            fu = 1.125_dp
        end if

        !! Set the He fudge factor
        !if ((Heswitch == 2) .or. (Heswitch == 5) .or. (Heswitch == 6)) then
        !    write(*,*)'Enter the fudge factor b_He'
        !    read(*,*)b_He
        !endif
        b_He = 0.86
    end subroutine set_switch
end module fudgefit


!##############################################################################
module cosmo_input
    use precision, only : dp
    use constants, only : pi, c, G, a_rad, kappa, parsec, m_1H, not4
    implicit none

    ! fnu is the contribution of neutrinos to the radn. energy density
    real(dp), parameter :: fnu = (21._dp / 8._dp) * (4._dp / 11._dp)**(4._dp / 3._dp)

    real(dp) :: H0      ! H0 is the Hubble constant in SI units [1/s]
    real(dp) :: OmegaM  ! Omega_0 contribution from total matter: OmegaM = OmegaC + OmegaB
    real(dp) :: OmegaL  ! Omega_0 contribution from a Cosmological constant
    real(dp) :: OmegaK  ! Omega_0 contribution in curvature (1 - O_T - O_L)
    real(dp) :: Tnow    ! Tnow is the observed CMB temperature today [K]
    real(dp) :: Nnow    ! Nnow is number density today
    real(dp) :: fHe     ! fHe is He/H number ratio = Yp / 4(1 - Yp)
    real(dp) :: z_eq    ! redshift of matter-radiation equality

    contains

    subroutine set_cosmo_input(OmegaB, OmegaC, OmegaL_in, H0_in, Tnow_in, Yp, Nnow_out, fHe_out)
        real(dp), intent(in)  :: OmegaB     ! Omega_0 in baryons today
        real(dp), intent(in)  :: OmegaC     ! Omega_0 in (cold) dark matter: OmegaM = OmegaC + OmegaB
        real(dp), intent(in)  :: OmegaL_in  ! input value for OmegaL
        real(dp), intent(in)  :: H0_in      ! input value of Hubble constant in units of 100 km/s/Mpc
        real(dp), intent(in)  :: Tnow_in    ! input value for Tnow
        real(dp), intent(in)  :: Yp         ! primordial helium abundace
        real(dp), intent(out) :: Nnow_out
        real(dp), intent(out) :: fHe_out

        real(dp) :: mu_H
        real(dp) :: mu_T

        ! convert the Hubble constant units
        H0 = H0_in * 1.e3_dp / (1.e6_dp * parsec)  ! from km/s/Mpc to 1/s
        Tnow = Tnow_in

        OmegaL = OmegaL_in
        OmegaM = OmegaB + OmegaC
        OmegaK = 1._dp - OmegaM - OmegaL  ! curvature

        ! sort out the helium abundance parameters
        ! mu_H, mu_T: mass per H atom and mass per particle
        mu_H = 1._dp / (1._dp - Yp)                 ! Mass per H atom
        mu_T = not4 / (not4 - (not4 - 1._dp) * Yp)  ! Mass per atom
        fHe = Yp / (not4 * (1._dp - Yp))            ! n_He_tot / n_H_tot
        fHe_out = fHe

        Nnow = 3._dp * H0**2 / kappa * OmegaB / (mu_H * m_1H)
        Nnow_out = Nnow

        ! (this is explictly for 3 massless neutrinos - change if N_nu /= 3)
        z_eq = 3._dp * H0**2 / kappa * c**2 / (a_rad * (1._dp + fnu) * Tnow**4) * (OmegaB + OmegaC)
        z_eq = z_eq - 1._dp
    end subroutine set_cosmo_input
end module cosmo_input


!##############################################################################
module recombination
    use precision, only : dp
    use constants, only : CB1_H, CB1_He1, CB1_He2, CR
    use fudgefit, only : set_switch
    use cosmo_input, only : set_cosmo_input
    use ode_solver, only: dverk
    implicit none
    contains

    subroutine recfast_func(OmegaB, OmegaC, OmegaL_in, H0_in, Tnow, Yp, Hswitch_in, Heswitch_in, &
                            zinitial, zfinal, tol, Nz, z_array, x_array)
        integer,  intent(in) :: Nz           ! number of output redshifts (integer)
        integer,  intent(in) :: Hswitch_in
        integer,  intent(in) :: Heswitch_in
        real(dp), intent(in) :: OmegaB
        real(dp), intent(in) :: OmegaC
        real(dp), intent(in) :: OmegaL_in
        real(dp), intent(in) :: H0_in
        real(dp), intent(in) :: Tnow
        real(dp), intent(in) :: Yp
        real(dp), intent(in) :: zinitial     ! starting redshift
        real(dp), intent(in) :: zfinal       ! ending redshift
        real(dp), intent(in) :: tol          ! tolerance for the integrator
        real(dp), intent(out) :: z_array(Nz)
        real(dp), intent(out) :: x_array(Nz)

        integer, parameter :: Ndim = 3  ! number of d.e.'s to solve (integer)

        integer :: i              ! loop index (integer)
        integer :: ind            ! ind, nw: work-space for dverk (integer)
        integer :: nw             ! ind, nw: work-space for dverk (integer)
        real(dp) :: Trad, Tmat
        real(dp) :: x_H0, x_He0
        real(dp) :: x0
        real(dp) :: w0, w1, logw0, logw1, hw
        real(dp) :: z_old, z_new  ! z_old and z_new are for each pass to the integrator
        real(dp) :: rhs           ! rhs is dummy for calculating x0
        real(dp) :: x_H, x_He

        real(dp) :: cw(24)      ! cw(24), w(3,9): work space for dverk
        real(dp) :: w(Ndim, 9)  ! cw(24), w(3,9): work space for dverk

        ! DLSODE, DLSODA, and VODE parameters
        integer  :: neq(1)
        integer  :: itol = 1
        real(dp) :: atol(1)
        integer  :: itask = 1
        integer  :: istate = 1
        integer  :: iopt = 0
        integer  :: jt = 2
        integer, parameter :: lrw = 22 + Ndim * max(16, Ndim + 9)  ! length of rwork
        integer, parameter :: liw = 20 + Ndim                      ! length of iwork
        real(dp) :: rwork(lrw)  ! work space for DLSODA
        integer  :: iwork(liw)  ! work space for DLSODA
        !real(dp) :: rpar
        !integer  :: ipar
        integer :: mf = 22
        real(dp) :: jdummy

        real(dp) :: y(Ndim)

        real(dp) :: Nnow
        real(dp) :: fHe

        external ion

        call set_cosmo_input(OmegaB, OmegaC, OmegaL_in, H0_in, Tnow, Yp, Nnow, fHe)

        call set_switch(Hswitch_in, Heswitch_in)

        ! Set initial matter temperature
        Tmat = Tnow * (1._dp + zinitial)  ! Initial rad. & mat. temperature

        call get_init(zinitial, x_H0, x_He0, x0)

        ! OK that's the initial conditions, now start computing our arrays
        w0 = 1._dp / sqrt(1._dp + zinitial)       ! conformal-time-like initial zi
        w1 = 1._dp / sqrt(1._dp + zfinal)         ! conformal-time-like final zf
        logw0 = log(w0)                           ! log of w0
        logw1 = log(w1)                           ! log of w1
        hw = (logw1 - logw0) / real(Nz, kind=dp)  ! interval in log of conf time

        ! Set up work-space stuff for dverk
        ind = 1
        nw  = Ndim
        do i = 1, 24
            cw(i) = 0._dp
        end do

        y(1) = x_H0   ! ionised fraction of H  - y(1) in R-K routine
        y(2) = x_He0  ! ionised fraction of He - y(2) in R-K routine, (note that x_He=n_He+/n_He here and not n_He+/n_H)
        y(3) = Tmat   ! matter temperature     - y(3) in R-K routine

        do i = 1, Nz
            ! calculate the start and end redshift for the interval at each z
            ! or just at each z
            z_old = zinitial + real(i-1, kind=dp) * (zfinal - zinitial) / real(Nz, kind=dp)
            z_new = zinitial + real(i  , kind=dp) * (zfinal - zinitial) / real(Nz, kind=dp)

            ! Use Saha to get x_e, using the equation for x_e for ionised helium
            ! and for neutral helium.
            ! Everyb_trip ionised above z=8000.  First ionization over by z=5000.
            ! Assume He all singly ionised down to z=3500, then use He Saha until
            ! He is 99% singly ionised, and *then* switch to joint H/He recombination.
            if (z_new > 8000._dp) then  ! everything ionised: x_H0=1, x_He0=1
                x_H0 = 1._dp
                x_He0 = 1._dp
                x0 = 1._dp + 2._dp * fHe
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Tnow * (1._dp + z_new)
            else if (z_new > 5000._dp) then
                x_H0 = 1._dp
                x_He0 = 1._dp
                rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z_new)) - CB1_He2 / (Tnow * (1._dp + z_new))) / Nnow
                rhs = rhs * 1._dp  ! ratio of g's is 1 for He++ <-> He+
                x0 = 0.5_dp * (sqrt((rhs - 1._dp - fHe)**2 + (4._dp + 8._dp * fHe) * rhs) - (rhs - 1._dp - fHe))
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Tnow * (1._dp + z_new)
            else if (z_new > 3500._dp) then
                x_H0 = 1._dp
                x_He0 = 1._dp
                x0 = x_H0 + fHe * x_He0
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Tnow * (1._dp + z_new)
            else if (y(2) > 0.99_dp) then
                x_H0 = 1._dp
                rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z_new)) - CB1_He1 / (Tnow * (1._dp + z_new))) / Nnow
                rhs = rhs * 4._dp      !ratio of g's is 4 for He+ <-> He0
                x_He0 = 0.5_dp * (sqrt( (rhs - 1._dp)**2 + 4._dp * (1._dp + fHe) * rhs ) - (rhs - 1._dp))
                x0 = x_He0
                x_He0 = (x0 - 1._dp) / fHe
                y(1) = x_H0
                y(2) = x_He0
                y(3) = Tnow * (1._dp + z_new)
!            else if (y(1) > 0.99_dp) then
!                write(*, '(i1.0,2x,f9.2,2x,f9.2,2x,e24.18,2x,e24.18,2x,e24.18)') 1, z_old, z_new, y(1), y(2), y(3)
!                rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z_new)) - CB1_H / (Tnow * (1._dp + z_new))) / Nnow
!                x_H0 = 0.5_dp * (sqrt( rhs**2 + 4._dp * rhs ) - rhs )
!!                call dverk(nw, ion, z_old, y, z_new, tol, ind, cw, nw, w)
!                call DLSODA(ion, Ndim, y, z_old, z_new, itol, tol, tol, itask, istate, iopt, rwork, lrw, iwork, liw, jdummy, jt)
!                y(1) = x_H0
!                x0 = y(1) + fHe * y(2)
!                write(*, '(i1.0,2x,f9.2,2x,f9.2,2x,e24.18,2x,e24.18,2x,e24.18)') 1, z_old, z_new, y(1), y(2), y(3)
!                write(*, *) istate
            else
                write(*, '(i1.0,2x,f9.2,2x,f9.2,2x,e24.18,2x,e24.18,2x,e24.18)') 2, z_old, z_new, y(1), y(2), y(3)
!                call dverk(nw, ion, z_old, y, z_new, tol, ind, cw, nw, w)
                call DLSODA(ion, Ndim, y, z_old, z_new, itol, tol, tol, itask, istate, iopt, rwork, lrw, iwork, liw, jdummy, jt)
                x0 = y(1) + fHe * y(2)
                write(*, '(i1.0,2x,f9.2,2x,f9.2,2x,e24.18,2x,e24.18,2x,e24.18)') 2, z_old, z_new, y(1), y(2), y(3)
                write(*, *) istate
            end if

            Trad = Tnow * (1._dp + z_new)  ! Trad and Tmat are radiation and matter temperatures
            x_H = y(1)   ! ionised fraction of H  - y(1) in R-K routine
            x_He = y(2)  ! ionised fraction of He - y(2) in R-K routine, (note that x_He=n_He+/n_He here and not n_He+/n_H)
            Tmat = y(3)  ! matter temperature     - y(3) in R-K routine

            z_array(i) = z_new
            x_array(i) = x0

        end do
    end subroutine recfast_func
end module recombination


!##############################################################################
program recfast
    use precision, only : dp
    use recombination, only : recfast_func
    implicit none

!   --- Arguments
    integer :: i  ! loop index (integer)
    integer :: Nz
    integer Heswitch_in, Hswitch_in

    real(dp) :: tol
    real(dp) :: zinitial, zfinal
    real(dp) :: OmegaB, OmegaC, OmegaL
    real(dp) :: H0_in, Tnow, Yp

    real(dp), dimension(:), allocatable :: z_array  ! array of redshifts written to file in the end
    real(dp), dimension(:), allocatable :: x_array  ! array of ionised fraction written to file in the end

    character(len=80) :: fileout
    
    tol = 1.e-8_dp  ! tolerance for R-K

    !   ###########################################################################
    write(*,*)'recfast version 1.8'
    write(*,*)'Using Hummer''s case B recombination rates for H'
    write(*,*)' with H fudge factor = 1.14 (or 1.125 plus high z fit),'
    write(*,*)' b_He fudge factor = 0.86,'
    write(*,*)' and a fit to tabulated HeII singlet recombination rates'
    write(*,*)

    ! will output every 10 in z, but this is easily changed also

    write(*,*) 'Enter output file name'
    read(*,'(a)') fileout
    write(*,*)
    write(*,*) 'Enter zinital, zfinal, Nz (e.g. 10000 0 1000)'
    read(*,*) zinitial, zfinal, Nz
    allocate(z_array(Nz))
    allocate(x_array(Nz))
    write(*,*)
    write(*,*) 'Enter Omega_B, Omega_DM, Omega_vac (e.g. 0.04 0.20 0.76)'
    read(*,*) OmegaB, OmegaC, OmegaL
    !write(*,'(1x,''Omega_K = '',f4.2)') OmegaK
    write(*,*)
    write(*,*) 'Enter H_0 (in km/s/Mpc), T_0, Y_p (e.g. 70 2.725 0.25)'
    read(*,*) H0_in, Tnow, Yp
    write(*,*)
    ! Modification for H correction (Hswitch):
    write(*,*) 'Modification for H recombination:'
    write(*,*) '0) no change from old Recfast'
    write(*,*) '1) include correction'
    write(*,*) 'Enter the choice of modification for H (0-1):'
    read(*,*) Hswitch_in
    write(*,*)
    ! Modification for HeI recombination (Heswitch):
    write(*,*) 'Modification for HeI recombination:'
    write(*,*) '0) no change from old Recfast'
    write(*,*) '1) full expression for escape probability for singlet'
    write(*,*) '   1P-1S transition'
    write(*,*) '2) also including effect of contiuum opacity of H on HeI'
    write(*,*) '   singlet (based in fitting formula suggested by'
    write(*,*) '   Kholupenko, Ivanchik & Varshalovich, 2007)'
    write(*,*) '3) only including recombination through the triplets'
    write(*,*) '4) including 3 and the effect of the contiuum '
    write(*,*) '   (although this is probably negligible)'
    write(*,*) '5) including only 1, 2 and 3'
    write(*,*) '6) including all of 1 to 4'
    write(*,*) 'Enter the choice of modification for HeI (0-6):'
    read(*,*) Heswitch_in
    write(*,*)

    ! OK that's the initial conditions, now start writing output file
    call recfast_func(OmegaB, OmegaC, OmegaL, H0_in, Tnow, Yp, Hswitch_in, Heswitch_in, &
                      zinitial, zfinal, tol, Nz, z_array, x_array)

    open(unit=7, status='new', form='formatted', file=fileout)
    write(7, '(''#'', 1x,''z'', 2x, ''x_e'')')
    do i = 1, Nz
        write(7, '(f9.2,2x,e24.18)') z_array(i), x_array(i)
    end do

end program recfast


!##############################################################################
subroutine get_init(z, x_H0, x_He0, x0)
!   Set up the initial conditions so it will work for general,
!   but not pathological choices of z_old
!   Initial ionization fraction using Saha for relevant species
    use precision, only : dp
    use constants, only : CB1_H, CB1_He1, CB1_He2, CR
    use cosmo_input, only : Tnow, Nnow, fHe
    implicit none

    real(dp), intent(in) :: z       ! initial redshift
    real(dp), intent(out) :: x_H0   ! initial ionised fraction of Hydrogen
    real(dp), intent(out) :: x_He0  ! initial ionised fraction of Helium
    real(dp), intent(out) :: x0     ! initial ionised fraction
    real(dp) :: rhs

    if (z > 8000._dp) then
        x_H0 = 1._dp
        x_He0 = 1._dp
        x0 = 1._dp + 2._dp * fHe
    else if (z > 3500._dp) then
        x_H0 = 1._dp
        x_He0 = 1._dp
        rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z)) - CB1_He2 / (Tnow * (1._dp + z))) / Nnow
        rhs = rhs * 1._dp      !ratio of g's is 1 for He++ <-> He+
        x0 = 0.5_dp * (sqrt((rhs - 1._dp - fHe)**2 + (4._dp + 8._dp * fHe) * rhs) - (rhs - 1._dp - fHe))
    else if (z > 2000._dp) then
        x_H0 = 1._dp
        rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z)) - CB1_He1 / (Tnow * (1._dp + z))) / Nnow
        rhs = rhs * 4._dp      !ratio of g's is 4 for He+ <-> He0
        x_He0 = 0.5_dp * (sqrt((rhs - 1._dp)**2 + 4._dp * (1._dp + fHe) * rhs) - (rhs - 1._dp))
        x0 = x_He0
        x_He0 = (x0 - 1._dp) / fHe
    else
        rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z)) - CB1_H / (Tnow * (1._dp + z))) / Nnow
        x_H0 = 0.5_dp * (sqrt( rhs**2 + 4._dp * rhs) - rhs)
        x_He0 = 0._dp
        x0 = x_H0
    end if
    return
end subroutine get_init


!##############################################################################
subroutine ion(Ndim, z, y, f)
    use precision, only : dp
    use constants, only : pi, c, h_P, k_B
    use constants, only : m_1H, not4, Lambda_H, Lambda_He, L_He_2p
    use constants, only : A2P_s, A2P_t, sigma_He_2Ps, sigma_He_2Pt, L_He_2Pt, L_He_2St, L_He2St_ion
    use constants, only : CDB_H, CDB_He, CR, CK_H, CK_He, CL_H, CL_He, CT, Bfact
    use fudgefit, only : AGauss1, AGauss2, zGauss1, zGauss2, wGauss1, wGauss2
    use fudgefit, only : H_frac, Hswitch, Heswitch, fu, b_He
    use cosmo_input, only : OmegaL, OmegaM, OmegaK, H0, Tnow, Nnow, fHe, z_eq
    implicit none

    integer,  intent(in)  :: Ndim      ! number of d.e.'s to solve (integer)
    real(dp), intent(in)  :: z         ! redshift - w is sqrt(1 + z), like conformal time
    real(dp), intent(in)  :: y(Ndim)   ! y = [x_H, x_He, Tmat]
    real(dp), intent(out) :: f(Ndim)   ! f's are the derivatives of the y's

    integer Heflag

    real(dp) :: x, nd_H, n_He, Trad, Tmat, x_H, x_He
    real(dp) :: Hz
    real(dp) :: Rup, Rdown, K, K_He, Rup_He, Rdown_He, He_Boltz
    real(dp) :: timeTh, timeH
    real(dp) :: a_VF, b_VF, T_0, T_1, sq_0, sq_1, a_PPB, b_PPB, c_PPB, d_PPB
    real(dp) :: tauHe_s, pHe_s
    real(dp) :: Doppler, gamma_2Ps, pb, qb, AHcon
    real(dp) :: a_trip, b_trip, Rdown_trip, Rup_trip
    real(dp) :: tauHe_t, pHe_t, CL_PSt, CfHe_t, gamma_2Pt
    real(dp) :: dHdz, epsilon

    ! the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
    a_PPB = 4.309_dp
    b_PPB = -0.6166_dp
    c_PPB = 0.6703_dp
    d_PPB = 0.5300_dp
    ! the Verner and Ferland type fitting parameters for Helium
    ! fixed to match those in the SSS papers, and now correct
    a_VF = 10._dp**(-16.744_dp)
    b_VF = 0.711_dp
    T_0 = 10._dp**(0.477121_dp)   !3K
    T_1 = 10._dp**(5.114_dp)
    ! fitting parameters for HeI triplets
    ! (matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)

    ! Parameters and quantities to describe the extra triplet states
    ! and also the continuum opacity of H, with a fitting function
    ! suggested by KIV, astro-ph/0703438
    a_trip = 10._dp**(-16.306_dp)  ! used to fit HeI triplet recombination rate
    b_trip = 0.761_dp              ! used to fit HeI triplet recombination rate

    x_H = y(1)            ! ionised fraction of H  - y(1) in R-K routine
    x_He = y(2)           ! ionised fraction of He - y(2) in R-K routine, (note that x_He=n_He+/n_He here and not n_He+/n_H)
    Tmat = y(3)           ! matter temperature     - y(3) in R-K routine
    x = x_H + fHe * x_He  ! total ionised fraction, relative to H

    nd_H = Nnow * (1._dp + z)**3
    n_He = fHe * Nnow * (1._dp + z)**3
    Trad = Tnow * (1._dp + z)
    ! Hz is the value of the Hubble parameter H at the specific z
    Hz = H0 * sqrt(OmegaM * (1._dp + z)**4 / (1._dp + z_eq) &
                   + OmegaM * (1._dp + z)**3 &
                   + OmegaK * (1._dp + z)**2 &
                   + OmegaL)

    ! dHdz is the derivative of the Hubble parameter H at the specific z
    dHdz = (H0**2 / 2._dp / Hz) * (  4._dp * OmegaM * (1._dp + z)**3 / (1._dp + z_eq) &
                                   + 3._dp * OmegaM * (1._dp + z)**2 &
                                   + 2._dp * OmegaK * (1._dp + z))

    ! Get the radiative rates using PPQ fit (identical to Hummer's table)
    Rdown = 1.e-19_dp * a_PPB * (Tmat / 1.e4_dp)**b_PPB / (1._dp + c_PPB * (Tmat / 1.e4_dp)**d_PPB)
    Rup = Rdown * (CR * Tmat)**(1.5_dp) * exp(-CDB_H / Tmat)

    ! calculate He using a fit to a Verner & Ferland type formula
    sq_0 = sqrt(Tmat / T_0)
    sq_1 = sqrt(Tmat / T_1)
    ! typo here corrected by Wayne Hu and Savita Gahlaut
    Rdown_He = a_VF / (sq_0 * (1._dp + sq_0)**(1._dp - b_VF))
    Rdown_He = Rdown_He / (1._dp + sq_1)**(1._dp + b_VF)
    Rup_He = Rdown_He * (CR * Tmat)**(1.5_dp) * exp(-CDB_He / Tmat)
    Rup_He = 4._dp * Rup_He    !statistical weights factor for HeI
    ! Avoid overflow (pointed out by Jacques Roland)
    if((Bfact / Tmat) > 680._dp)then
        He_Boltz = exp(680._dp)
    else
        He_Boltz = exp(Bfact / Tmat)
    end if

    ! now deal with H and its fudges
    if (Hswitch == 0) then
        ! use Peebles coefficient K=lambda_a^3/8piH
        K = CK_H / Hz
    else
        ! fit a double Gaussian correction function
        K = CK_H / Hz * (1.0_dp + AGauss1 * exp(-((log(1.0_dp + z) - zGauss1) / wGauss1)**2) &
                                + AGauss2 * exp(-((log(1.0_dp + z) - zGauss2) / wGauss2)**2))
    end if

    ! add the HeI part, using same T_0 and T_1 values
    Rdown_trip = a_trip / (sq_0 * (1._dp + sq_0)**(1.0 - b_trip))
    Rdown_trip = Rdown_trip / ((1._dp + sq_1)**(1._dp + b_trip))
    Rup_trip = Rdown_trip * exp(-h_P * c * L_He2St_ion / (k_B * Tmat))
    Rup_trip = Rup_trip * ((CR * Tmat)**(1.5_dp)) * (4._dp / 3._dp)
    ! last factor here is the statistical weight

    ! try to avoid "NaN" when x_He gets too small
    if ((x_He < 5.e-9_dp) .or. (x_He > 0.980)) then
        Heflag = 0
    else
        Heflag = Heswitch
    end if
    if (Heflag == 0) then
        ! use Peebles coefficient for He
        K_He = CK_He / Hz
    else  ! for Heflag>0       
        ! use Sobolev escape probability
        tauHe_s = A2P_s * CK_He * 3._dp * n_He * (1._dp - x_He) / Hz
        pHe_s = (1._dp - exp(-tauHe_s)) / tauHe_s
        K_He = 1._dp / (A2P_s * pHe_s * 3._dp * n_He * (1._dp - x_He))
        if (((Heflag == 2) .or. (Heflag >= 5)) .and. (x_H < 0.9999999_dp)) then
            ! smoother criterion here from Antony Lewis & Chad Fendt
            ! use fitting formula for continuum opacity of H
            ! first get the Doppler width parameter
            Doppler = 2._dp * k_B * Tmat / (m_1H * not4 * c * c)
            Doppler = c * L_He_2p * sqrt(Doppler)
            gamma_2Ps = 3._dp * A2P_s * fHe * (1._dp - x_He) * c * c &
                        / (sqrt(pi) * sigma_He_2Ps * 8._dp * pi * Doppler * (1._dp - x_H)) / ((c * L_He_2p)**2)
            pb = 0.36_dp  ! value from KIV (2007)
            qb = b_He
            ! calculate AHcon, the value of A * p_(con, H) for H continuum opacity
            AHcon = A2P_s / (1._dp + pb * (gamma_2Ps**qb))
            K_He = 1._dp / ((A2P_s * pHe_s + AHcon) * 3._dp * n_He * (1._dp - x_He))
        end if
        if (Heflag >= 3) then     !include triplet effects
            tauHe_t = A2P_t * n_He * (1._dp - x_He) * 3._dp
            tauHe_t = tauHe_t /(8._dp * pi * Hz * L_He_2Pt**3)
            pHe_t = (1._dp - exp(-tauHe_t)) / tauHe_t
            CL_PSt = h_P * c * (L_He_2Pt - L_He_2St) / k_B
            if ((Heflag == 3) .or. (Heflag == 5) .or. (x_H > 0.99999_dp)) then
                ! no H cont. effect
                ! CfHe_t: triplet statistical correction
                CfHe_t = A2P_t * pHe_t * exp(-CL_PSt / Tmat)
                CfHe_t = CfHe_t / (Rup_trip + CfHe_t)   !"C" factor for triplets
            else                  !include H cont. effect
                Doppler = 2._dp * k_B * Tmat / (m_1H * not4 * c * c)
                Doppler = c * L_He_2Pt * sqrt(Doppler)
                gamma_2Pt = 3._dp * A2P_t * fHe * (1._dp - x_He) * c * c &
                            / (sqrt(pi) * sigma_He_2Pt * 8._dp * pi * Doppler * (1._dp - x_H) * (c * L_He_2Pt)**2)
                ! use the fitting parameters from KIV (2007) in this case
                pb = 0.66_dp
                qb = 0.9_dp
                AHcon = A2P_t / (1._dp + pb * gamma_2Pt**qb) / 3._dp
                CfHe_t = (A2P_t * pHe_t + AHcon) * exp(-CL_PSt / Tmat)
                CfHe_t = CfHe_t / (Rup_trip + CfHe_t)   !"C" factor for triplets
            end if
        end if
    end if

    ! Estimates of Thomson scattering time and Hubble time
    timeTh = (1._dp / (CT * Trad**4)) * (1._dp + x + fHe) / x   !Thomson time
    timeH = 2._dp / (3._dp * H0 * (1._dp + z)**(1.5_dp))        !Hubble time

    ! calculate the derivatives
    ! turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    ! (clunky, but seems to work)
    if (x_H > 0.99_dp) then         !don't change at all
        f(1) = 0._dp
    !else if ((x_H > 0.98_dp) .and. (Heflag == 0)) then    !don't modify
    else if (x_H > 0.985_dp) then     !use Saha rate for Hydrogen
        f(1) = (x * x_H * nd_H * Rdown - Rup * (1._dp - x_H) * exp(-CL_H / Tmat)) / (Hz * (1._dp + z))
        ! AL: following commented as not used
        !! for interest, calculate the correction factor compared to Saha
        !! (without the fudge)
        !factor = (1._dp + K * Lambda_H * nd_H * (1._dp - x_H)) &
        !         / (Hz * (1._dp + z) * (1._dp + K * Lambda_H * nd_H * (1._dp - x) + K * Rup * nd_H * (1._dp - x)))
    else  ! use full rate for H
        f(1) = ((x * x_H * nd_H * Rdown - Rup * (1._dp - x_H) * exp(-CL_H / Tmat)) &
                * (1._dp + K * Lambda_H * nd_H * (1._dp - x_H))) &
               / (Hz * (1._dp+z) * (1._dp/fu + K * Lambda_H * nd_H * (1._dp-x_H) / fu + K * Rup * nd_H * (1._dp-x_H)))
    end if
    ! turn off the He once it is small
    if (x_He < 1.e-15_dp) then
        f(2) = 0._dp
    else
        f(2) = ((x * x_He * nd_H * Rdown_He - Rup_He * (1._dp - x_He) * exp(-CL_He / Tmat)) &
                * (1._dp + K_He * Lambda_He * n_He * (1._dp - x_He) * He_Boltz)) &
               / (Hz * (1._dp + z) * (1._dp + K_He * (Lambda_He + Rup_He) * n_He * (1._dp - x_He) * He_Boltz))
        ! Modification to HeI recombination including channel via triplets
        if (Heflag >= 3) then
            f(2) = f(2) + (x * x_He * nd_H * Rdown_trip - Rup_trip * (1._dp-x_He) * 3._dp * exp(-h_P * c * L_He_2St &
                                                                                                / (k_B * Tmat))) &
                          * CfHe_t / (Hz * (1._dp + z))
        end if
    end if

    ! follow the matter temperature once it has a chance of diverging

    if (timeTh < H_frac * timeH) then
        !f(3) = Tmat / (1._dp + z)  !Tmat follows Trad
        ! additional term to smooth transition to Tmat evolution,
        ! (suggested by Adam Moss)
        epsilon = Hz * (1._dp + x + fHe) / (CT * Trad**3 * x)  ! approximate difference (=Trad-Tmat) at high z
        f(3) = Tnow &
               + epsilon * ((1._dp + fHe) / (1._dp + fHe + x)) * ((f(1) + fHe * f(2)) / x) &
               - epsilon * dHdz / Hz + 3.0_dp * epsilon / (1._dp + z)
    else
        f(3) = CT * Trad**4 * x / (1._dp + x + fHe) * (Tmat - Trad) / (Hz * (1._dp + z)) + 2._dp * Tmat / (1._dp + z)
    end if
    return
end subroutine ion


!##############################################################################
! jacobian for DLSODA:
!subroutine jdummy(Ndim, z, y, ML, MU, pd, nrowpd)
!    use precision, only : dp
!    implicit none
!    integer  :: Ndim              ! number of d.e.'s to solve (integer)
!    integer  :: ML
!    integer  :: MU
!    integer  :: nrowpd
!    real(dp) :: z                 ! redshift - w is sqrt(1 + z), like conformal time
!    real(dp) :: y(Ndim)           ! y = [x_H, x_He, Tmat]
!    real(dp) :: pd(nrowpd, Ndim)  ! f's are the partial derivatives of the y's
!end subroutine jdummy
!
!SUBROUTINE FEX (NEQ, T, Y, YDOT)
!    use precision, only : dp
!    implicit none
!    integer  :: NEQ
!    real(dp) :: T
!    real(dp) :: Y(3)
!    real(dp) :: YDOT(3)
!    RETURN
!END SUBROUTINE FEX


!! example function for VODE:
!subroutine FEX (Ndim, z, y, f, rpar, ipar)
!    use precision, only : dp
!    implicit none
!    integer  :: Ndim      ! number of d.e.'s to solve (integer)
!    real(dp) :: z         ! redshift - w is sqrt(1 + z), like conformal time
!    real(dp) :: y(Ndim)   ! y = [x_H, x_He, Tmat]
!    real(dp) :: f(Ndim)   ! f's are the derivatives of the y's
!    real(dp) :: rpar
!    integer  :: ipar
!!    f(1) = -.04D0*y(1) + 1.D4*y(2)*y(3)
!!    f(3) = 3.D7*y(2)*y(2)
!!    f(2) = -f(1) - f(3)
!    return
!end


!! jacobian for VODE:
!subroutine JVODE (Ndim, z, y, ml, mu, pd, nrpd, rpar, ipar)
!    use precision, only : dp
!    implicit none
!    integer  :: Ndim      ! number of d.e.'s to solve (integer)
!    real(dp) :: z         ! redshift - w is sqrt(1 + z), like conformal time
!    real(dp) :: y(Ndim)   ! y = [x_H, x_He, Tmat]
!    integer  :: ml
!    integer  :: mu
!    integer  :: nrpd
!    real(dp) :: pd(nrpd, Ndim)
!    real(dp) :: rpar
!    integer  :: ipar
!    !    pd(1,1) = -.04D0
!    !    pd(1,2) = 1.D4*y(3)
!    !    pd(1,3) = 1.D4*y(2)
!    !    pd(2,1) = .04D0
!    !    pd(2,3) = -pd(1,3)
!    !    pd(3,2) = 6.D7*y(2)
!    !    pd(2,2) = -pd(1,2) - pd(3,2)
!    return
!end
