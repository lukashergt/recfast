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
!V  Version: 1.6.0
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
!D      imput their own H(z) to account for extra physics.
!D      This is provided as a PROGRAM, which can be easily converted
!D      to a SUBROUTINE for use in CMB Boltzmann codes.
!
!A  Arguments:
!A      Name, Description
!A      Double precision throughout
!A
!A      z is redshift - W is sqrt(1 + z), like conformal time
!A      x is total ionised fraction, relative to H
!A      x_H is ionized fraction of H - y(1) in R-K routine
!A      x_He is ionized fraction of He - y(2) in R-K routine
!A      (note that x_He=n_He+/n_He here and not n_He+/n_H)
!A      Tmat is matter temperature - y(3) in R-K routine
!A      f's are the derivatives of the Y's
!A      alphaB is case B recombination rate
!A      alpHe is the singlet only HeII recombination rate
!A      a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!A      b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!A      c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!A      d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
!A      a_VF is Verner and Ferland type fitting parameter for Helium
!A      b_VF is Verner and Ferland type fitting parameter for Helium
!A      T_0 is Verner and Ferland type fitting parameter for Helium
!A      T_1 is Verner and Ferland type fitting parameter for Helium
!A      Tnow is the observed CMB temperature today
!A      OmegaT is the total Omega_0
!A      OmegaL is the Omega_0 contribution from a Cosmological constant
!A      OmegaK is the Omega_0 contribution in curvature (1 - O_T - O_L)
!A      OmegaB is Omega in baryons today
!A      OmegaC is the Omega_0 in (cold) dark matter: OmegaT = OmegaC + OmegaB
!A      Yp is the primordial helium abundace
!A      fHe is He/H number ratio = Yp / 4(1 - Yp)
!A      Trad and Tmat are radiation and matter temperatures
!A      epsilon is the approximate difference (=Trad-Tmat) at high z
!A      OmegaB is Omega in baryons today
!A      H0inp is input value of Hubble constant in units of 100 km/s/Mpc
!A      H0 is Hubble constant in SI units
!A      Hz is the value of H at the specific z (in ion)
!A      G is grvitational constant
!A      n is number density of hydrogen
!A      Nnow is number density today
!A      x0 is initial ionized fraction
!A      x_H0 is initial ionized fraction of Hydrogen
!A      x_He0 is initial ionized fraction of Helium
!A      rhs is dummy for calculating x0
!A      zinitial and zfinal are starting and ending redshifts
!A      fnu is the contribution of neutrinos to the radn. energy density
!A      zeq is the redshift of matter-radiation equality
!A      zstart and zend are for each pass to the integrator
!A      w0 and w1 are conformal-time-like initial and final zi and zf's
!A      Lw0 and Lw1 are logs of w0 and w1
!A      hw is the interval in W
!A      mu_H, mu_T: mass per H atom and mass per particle
!A      H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
!A      dHdz is the derivative of H at the specific z (in ion)
!A      fu is a "fudge factor" for H, to approximate low z behaviour
!A      b_He is a "fudge factor" for HeI, to approximate higher z behaviour
!A      Heswitch is an integer for modifying HeI recombination
!A      Parameters and quantities to describe the extra triplet states
!A       and also the continuum opacity of H, with a fitting function
!A       suggested by KIV, astro-ph/0703438
!A      a_trip: used to fit HeI triplet recombination rate
!A      b_trip: used to fit HeI triplet recombination rate
!A      CL_PSt = h_P * c * (L_He_2Pt - L_He_2St) / k_B
!A      CfHe_t: triplet statistical correction
!A      Hswitch is an integer for modifying the H recombination
!A      AGauss1 is the amplitude of the 1st Gaussian for the H fudging
!A      AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
!A      zGauss1 is the ln(1 + z) central value of the 1st Gaussian
!A      zGauss2 is the ln(1 + z) central value of the 2nd Gaussian
!A      wGauss1 is the width of the 1st Gaussian
!A      wGauss2 is the width of the 2nd Gaussian
!A      tol: tolerance for the integrator
!A      cw(24), w(3,9): work space for dverk
!A      Ndim: number of d.e.'s to solve (integer)
!A      Nz: number of output redshitf (integer)
!A      I: loop index (integer)
!A      ind, nw: work-space for dverk (integer)
!
!G  Global data (common blocks) referenced:
!G      /zLIST/zinitial, zfinal, Nz
!G      /data/H_frac, fHe
!G      /Cosmo/Tnow, H0, Nnow, z_eq, OmegaT, OmegaL, OmegaK
!G      /Hemod/b_He
!G      /Hmod/AGauss1, AGauss2, zGauss1, zGauss2, wGauss1, wGauss2
!G      /Switch/Heswitch, Hswitch
!
!F  File & device access:
!F      Unit    /I, IO, O /Name (if known)
!
!M  Modules used:
!M      precision (declares `dp` for double precision)
!M      constants (module for mathematical and physical constants and chemical data)
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
!H      Jul 1998 (include OmegaT correction and H fudge factor)
!H      Nov 1998 (change Trad to Tmat in Rup)
!H      Jan 1999 (tidied up for public consumption)
!H      Sep 1999 (switch to formula for alpha's, fix glitch)
!H      Feb 2000 (fixed overflow problem in He_Boltz)
!H      Oct 2001 (fixed OmegaT in z_eq)
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
!
!##############################################################################


!##############################################################################
module precision
    use, intrinsic :: iso_fortran_env
    implicit none

    integer, parameter :: sp = REAL32
    integer, parameter :: dp = REAL64
    integer, parameter :: qp = REAL128
end module precision


!##############################################################################
module constants
    use precision, only : dp
    implicit none

    ! mathematical constants
    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp ! pi, Archimedes' constant

    ! Codata internationally recommended 2018 values (in SI units):
    real(dp), parameter :: c = 2.99792458e8_dp           ! speed of light in vacuum [m/s]
    real(dp), parameter :: G = 6.67430e-11_dp            ! Newton's constant of gravitation [m^3/kg/s^2]
    real(dp), parameter :: h_P = 6.62607015e-34_dp       ! Planck's constant [kg m^2/s]
    real(dp), parameter :: k_B = 1.380649e-23_dp         ! Boltzmann's constant [kg m^2/s^2/K]
    real(dp), parameter :: m_e = 9.1093837015e-31_dp     ! electron mass [kg]
    real(dp), parameter :: amu = 1.66053906660e-27_dp    ! atomic mass unit [kg]
    real(dp), parameter :: sigma_e = 6.6524587321e-29_dp ! Thomson cross section [m^2]

    ! black body radiation constant for u=aT^4:
    real(dp), parameter :: a = 8._dp * pi**5 * k_B**4 / (15._dp * c**3 * h_P**3) ! [J/m^3/K^4]

    ! some astronomical units:
    real(dp), parameter :: AU = 149597870700._dp          ! astronomical unit [m]; exact definition of the IAU
    real(dp), parameter :: parsec = 648000._dp / pi * AU  ! parsec [m]; exact definition of the IAU

    ! Atomic mass evaluation (AME) 2020:
    real(dp), parameter :: m_1H_u  = 1.007825031898_dp   ! atomic mass of 1H in atomic mass units [u]
    real(dp), parameter :: m_2H_u  = 2.014101777844_dp   ! atomic mass of 2H in atomic mass units [u]
    real(dp), parameter :: m_3H_u  = 3.016049281320_dp   ! atomic mass of 3H in atomic mass units [u]
    real(dp), parameter :: m_3He_u = 3.016029321967_dp   ! atomic mass of 3He in atomic mass units [u]
    real(dp), parameter :: m_4He_u = 4.002603254130_dp   ! atomic mass of 4He in atomic mass units [u]
    real(dp), parameter :: m_1H    = m_1H_u * amu        ! atomic mass of 1H [kg]
    real(dp), parameter :: m_4He   = m_4He_u * amu       ! atomic mass of 4He [kg]
    real(dp), parameter :: not4    = m_4He_u / m_1H_u    ! atomic mass ratio of 4He to 1H
    ! ("not4" pointed out by Gary Steigman)

    ! 2-photon rates in SI units:
    real(dp), parameter :: Lambda_H  =  8.2290619_dp       ! 2s-1s two-photon decay rate for Hydrogen [1/s]; Sommerfeldt et al (2020)
    real(dp), parameter :: Lambda_He = 50.93_dp            ! 2s-1s two photon decay rate for Helium [1/s]; Bondy, Morton & Drake (2020)
    ! Atomic levels in SI units (from NIST 2020):
    real(dp), parameter :: L_H_ion   = 1.0967877174307e7_dp  ! Hydrogen ionization level [1/m]; Kramida, Ralchenko, Reader, and NIST ASD (2020)
    real(dp), parameter :: L_H_alpha = 8.2259163e6_dp        ! Hydrogen Lyman alpha wavenumber [1/m]; Kramida (2010)
    real(dp), parameter :: L_He1_ion = 1.9831066637e7_dp     ! Helium I ionization level [1/m]; Kandula, Gohle, Pinkert, Ubachs, Eikema (2010)
    real(dp), parameter :: L_He2_ion = 4.389088785e7_dp      ! Helium II ionization level [1/m]; Johnson & Soff (1985), rescaled by NIST
    real(dp), parameter :: L_He_2s   = 1.66277440141e7_dp    ! Helium I 2s level [1/m]; Morton, Wu & Drake, CJP (2006)
    real(dp), parameter :: L_He_2p   = 1.71134896946e7_dp    ! Helium I 2p (21P1-11S0) level [1/m]; Morton, Wu & Drake, CJP (2006)
    ! Atomic data for HeI:
    real(dp), parameter :: A2P_s        = 1.798287e9_dp      ! Einstein A coefficient for He 21P1-11S0; Morton, Wu & Drake (2006)
    real(dp), parameter :: A2P_t        = 177.58_dp          ! Einstein A coefficient for He 23P1-11S0; Lach & Pachuski (2001)
    real(dp), parameter :: L_He_2Pt     = 1.690871466e7_dp   ! level for 23P012-11S0 [1/m], Drake & Morton (2007)
    real(dp), parameter :: L_He_2St     = 1.5985597526e7_dp  ! level for 23S1-11S0 [1/m], Drake & Morton (2007)
    real(dp), parameter :: L_He2St_ion  = 3.8454693845e6_dp  ! level for 23S1-continuum [1/m], Drake & Morton (2007)
    real(dp), parameter :: sigma_He_2Ps = 1.436289e-22_dp    ! H ionization x-section at HeI 21P1-11S0 freq. [m^2]; Hummer & Storey (1998)
    real(dp), parameter :: sigma_He_2Pt = 1.484872e-22_dp    ! H ionization x-section at HeI 23P1-11S0 freq. [m^2]; Hummer & Storey (1998)

    ! some derived constants:
    real(dp), parameter :: Lalpha = 1._dp / L_H_alpha                   ! Hydrogen Lyman alpha wavelength [m]
    real(dp), parameter :: Lalpha_He = 1._dp / L_He_2p                  ! Helium I 2p-1s wavelength [m]
    real(dp), parameter :: DeltaB_H = h_P * c * (L_H_ion - L_H_alpha)   ! energy of first excited state from continuum = 3.4eV
    real(dp), parameter :: CDB_H = DeltaB_H / k_B                       ! Constants derived from B1, B2, R
    real(dp), parameter :: DeltaB_He = h_P * c * (L_He1_ion - L_He_2s)  ! 2s, not 2p; energy of first excited state from cont. for He = 3.4eV
    real(dp), parameter :: CDB_He = DeltaB_He / k_B                     ! n=2-infinity for He in Kelvin
    real(dp), parameter :: CB1_H = h_P * c * L_H_ion / k_B              ! CB1 = CDB * 4.; Lalpha and sigma_Th, calculated
    real(dp), parameter :: CB1_He1 = h_P * c * L_He1_ion / k_B          ! CB1 for HeI ionization potential
    real(dp), parameter :: CB1_He2 = h_P * c * L_He2_ion / k_B          ! CB1 for HeII ionization potential
    real(dp), parameter :: CR = 2._dp * pi * (m_e / h_P) * (k_B / h_P)  !
    real(dp), parameter :: CK_H = Lalpha**3 / (8._dp * pi)              !
    real(dp), parameter :: CK_He = Lalpha_He**3 / (8._dp * pi)          !
    real(dp), parameter :: CL_H = c * h_P / (k_B * Lalpha)              !
    real(dp), parameter :: CL_He = c * h_P / (k_B / L_He_2s)            ! comes from det.bal. of 2s-1s
    real(dp), parameter :: CT = (8._dp / 3._dp) * a * sigma_e / m_e / c !
    real(dp), parameter :: Bfact = h_P * c * (L_He_2p - L_He_2s) / k_B  ! Bfact = exp((E_2p - E_2s) / kT); Extra Boltzmann factor
end module constants


!##############################################################################
program recfast
    use precision, only : dp
    use constants, only : pi, c, G, a, parsec
    use constants, only : m_1H, not4
    use constants, only : CB1_H, CB1_He1, CB1_He2, CR
    implicit none

!   --- Arguments
    real(dp) :: Trad, Tmat
    real(dp) :: OmegaT, OmegaB, H0, H0inp, OmegaL, OmegaK, OmegaC
    real(dp) :: z, n, x, x0, rhs, x_H, x_He, x_H0, x_He0
    real(dp) :: Tnow, zinitial, zfinal, Nnow, z_eq, fnu
    real(dp) :: zstart, zend, w0, w1, Lw0, Lw1, hw
    real(dp) :: mu_H, mu_T, H_frac
    real(dp) :: Yp, fHe, fu, b_He
    real(dp) :: AGauss1, AGauss2, zGauss1, zGauss2, wGauss1, wGauss2

    real(dp) :: tol
    real(dp) :: cw(24), w(3,9)
    real(dp) :: y(3)

    integer Ndim, Nz, I
    integer ind, nw
    integer Heswitch, Hswitch

    character(len=80) :: fileout

!   --- Parameter statements
    parameter(tol = 1.e-5_dp)                !Tolerance for R-K

    external ion

!   --- Commons
    common/zLIST/zinitial, zfinal, Nz
    common/Cdata/H_frac, fHe, fu
    common/Hemod/b_He
    common/Hmod/AGauss1, AGauss2, zGauss1, zGauss2, wGauss1, wGauss2
    common/Switch/Heswitch, Hswitch

    common/Cosmo/Tnow, H0, Nnow, z_eq, OmegaT, OmegaL, OmegaK
!   ###########################################################################

!   --- Data
    data    AGauss1     /-0.14_dp/   !Amplitude of 1st Gaussian
    data    AGauss2     /0.079_dp/   !Amplitude of 2nd Gaussian
    data    zGauss1     /7.28_dp/    !ln(1 + z) of 1st Gaussian
    data    zGauss2     /6.73_dp/    !ln(1 + z) of 2nd Gaussian
    data    wGauss1     /0.18_dp/    !Width of 1st Gaussian
    data    wGauss2     /0.33_dp/    !Width of 2nd Gaussian
    ! Gaussian fits for extra H physics (fit by Adam Moss, modified by
    ! Antony Lewis)

    ! dimensions for integrator
    Ndim = 3

    write(*,*)'recfast version 1.5'
    write(*,*)'Using Hummer''s case B recombination rates for H'
    write(*,*)' with H fudge factor = 1.14 (or 1.125 plus high z fit),'
    write(*,*)' b_He fudge factor = 0.86,'
    write(*,*)' and a fit to tabulated HeII singlet recombination rates'
    write(*,*)

    ! These are easy to inquire as input, but let's use simple values
    zinitial = 1.e4_dp
    z = zinitial
    zfinal = 0._dp
    ! will output every 10 in z, but this is easily changed also

    write(*,*)'Enter output file name'
    read(*,'(a)')fileout

    write(*,*)'Enter Omega_B, Omega_DM, Omega_vac (e.g. 0.04 0.20 0.76)'
    read(*,*)OmegaB, OmegaC, OmegaL
    OmegaT = OmegaC + OmegaB            !total dark matter + baryons
    OmegaK = 1._dp - OmegaT - OmegaL   !curvature
    write(*,'(1x,''Omega_K = '',f4.2)')OmegaK
    write(*,*)
    write(*,*)'Enter H_0 (in km/s/Mpc), T_0, Y_p (e.g. 70 2.725 0.25)'
    read(*,*)H0inp, Tnow, Yp

    ! convert the Hubble constant units
    H0 = H0inp * 1.e3_dp / (1.e6_dp * parsec)  ! from km/s/Mpc to 1/s

    ! sort out the helium abundance parameters
    mu_H = 1._dp / (1._dp - Yp)           !Mass per H atom
    mu_T = not4 / (not4 - (not4 - 1._dp) * Yp)   !Mass per atom
    fHe = Yp / (not4 * (1._dp - Yp))       !n_He_tot / n_H_tot

    Nnow = 3._dp * H0 * H0 * OmegaB / (8._dp * pi * G * mu_H * m_1H)
    n = Nnow * (1._dp + z)**3
    fnu = (21._dp / 8._dp) * (4._dp / 11._dp)**(4._dp / 3._dp)
    ! (this is explictly for 3 massless neutrinos - change if N_nu /= 3)
    z_eq = (3._dp * (H0 * c)**2 / (8._dp * pi * G * a * (1._dp + fnu) * Tnow**4)) * OmegaT
    z_eq = z_eq - 1._dp

    ! Matter departs from radiation when t(Th) > H_frac * t(H)
    ! choose some safely small number
    H_frac = 1.e-3_dp

    ! Modification for H correction (Hswitch):
    write(*,*) 'Modification for H recombination:'
    write(*,*)'0) no change from old Recfast'
    write(*,*)'1) include correction'
    write(*,*)'Enter the choice of modification for H (0-1):'
    read(*,*)Hswitch

    ! Fudge factor to approximate the low z out of equilibrium effect
    if (Hswitch == 0) then
        fu = 1.14_dp
    else
        fu = 1.125_dp
    end if

    ! Modification for HeI recombination (Heswitch):
    write(*,*)'Modification for HeI recombination:'
    write(*,*)'0) no change from old Recfast'
    write(*,*)'1) full expression for escape probability for singlet'
    write(*,*)'   1P-1S transition'
    write(*,*)'2) also including effect of contiuum opacity of H on HeI'
    write(*,*)'   singlet (based in fitting formula suggested by'
    write(*,*)'   Kholupenko, Ivanchik & Varshalovich, 2007)'
    write(*,*)'3) only including recombination through the triplets'
    write(*,*)'4) including 3 and the effect of the contiuum '
    write(*,*)'   (although this is probably negligible)'
    write(*,*)'5) including only 1, 2 and 3'
    write(*,*)'6) including all of 1 to 4'
    write(*,*)'Enter the choice of modification for HeI (0-6):'
    read(*,*)Heswitch

!   Set the He fudge factor
!c  if((Heswitch == 2) .or. (Heswitch == 5) .or. (Heswitch == 6))then
!c    write(*,*)'Enter the fudge factor b_He'
!c    read(*,*)b_He
!c  endif
    b_He = 0.86

    ! Set initial matter temperature
    y(3) = Tnow * (1._dp + z)            !Initial rad. & mat. temperature
    Tmat = y(3)

    call get_init(z, x_H0, x_He0, x0)

    y(1) = x_H0
    y(2) = x_He0

    ! OK that's the initial conditions, now start writing output file

    open(unit=7, status='new', form='formatted', file=fileout)
    write(7, '(1x,''  z    '', 1x, ''     x_e   '')')

    w0 = 1._dp / sqrt(1._dp + zinitial) !like a conformal time
    w1 = 1._dp / sqrt(1._dp + zfinal)
    Lw0 = log(w0)
    Lw1 = log(w1)
    Nz = 1000
    hW = (Lw1 - Lw0) / real(Nz, kind=dp)     !interval in log of conf time

    ! Set up work-space stuff for dverk
    ind  = 1
    nw   = 3
    do i = 1, 24
        cw(i) = 0._dp
    end do

    do i = 1, Nz
        ! calculate the start and end redshift for the interval at each z
        ! or just at each z
        zstart = zinitial + real(i - 1, kind=dp) * (zfinal - zinitial) / real(Nz, kind=dp)
        zend   = zinitial + real(i, kind=dp) * (zfinal - zinitial) / real(Nz, kind=dp)

        ! Use Saha to get x_e, using the equation for x_e for ionized helium
        ! and for neutral helium.
        ! Everyb_trip ionized above z=8000.  First ionization over by z=5000.
        ! Assume He all singly ionized down to z=3500, then use He Saha until
        ! He is 99% singly ionized, and *then* switch to joint H/He recombination.

        z = zend

        if (zend > 8000._dp) then

            x_H0 = 1._dp
            x_He0 = 1._dp
            x0 = 1._dp + 2._dp * fHe
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow * (1._dp + z)

        else if(z > 5000._dp)then

            x_H0 = 1._dp
            x_He0 = 1._dp
            rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z)) - CB1_He2 / (Tnow * (1._dp + z))) / Nnow
            rhs = rhs * 1._dp      !ratio of g's is 1 for He++ <-> He+
            x0 = 0.5_dp * (sqrt((rhs - 1._dp - fHe)**2 + (4._dp + 8._dp * fHe) * rhs) - (rhs - 1._dp - fHe))
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow * (1._dp + z)

        else if(z > 3500._dp)then

            x_H0 = 1._dp
            x_He0 = 1._dp
            x0 = x_H0 + fHe * x_He0
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow * (1._dp + z)

        else if(y(2) > 0.99)then

            x_H0 = 1._dp
            rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z)) - CB1_He1 / (Tnow * (1._dp + z))) / Nnow
            rhs = rhs * 4._dp      !ratio of g's is 4 for He+ <-> He0
            x_He0 = 0.5_dp * (sqrt( (rhs - 1._dp)**2 + 4._dp * (1._dp + fHe) * rhs ) - (rhs - 1._dp))
            x0 = x_He0
            x_He0 = (x0 - 1._dp) / fHe
            y(1) = x_H0
            y(2) = x_He0
            y(3) = Tnow * (1._dp + z)

        else if (y(1) > 0.99_dp) then

            rhs = exp(1.5_dp * log(CR * Tnow / (1._dp + z)) - CB1_H / (Tnow * (1._dp + z))) / Nnow
            x_H0 = 0.5_dp * (sqrt( rhs**2 + 4._dp * rhs ) - rhs )

            call dverk(nw, ion, zstart, y, zend, tol, ind, cw, nw, w)
            y(1) = x_H0
            x0 = y(1) + fHe * y(2)

        else

            call dverk(nw, ion, zstart, y, zend, tol, ind, cw, nw, w)

            x0 = y(1) + fHe * y(2)

        end if

        Trad = Tnow * (1._dp + zend)
        Tmat = y(3)
        x_H = y(1)
        x_He = y(2)
        x = x0

        write(7, '(1x,f8.2,2x,g15.8)') zend, x

    end do

end program recfast


!##############################################################################
subroutine get_init(z, x_H0, x_He0, x0)
!   Set up the initial conditions so it will work for general,
!   but not pathological choices of zstart
!   Initial ionization fraction using Saha for relevant species
    use precision, only : dp
    use constants, only : CB1_H, CB1_He1, CB1_He2, CR
    implicit none

    real(dp) :: OmegaT, H0, OmegaL, OmegaK
    real(dp) :: z, x0, rhs, x_H0, x_He0
    real(dp) :: Tnow, Nnow, z_eq
    real(dp) :: H_frac
    real(dp) :: fHe, fu

    common/Cdata/H_frac, fHe, fu
    common/Cosmo/Tnow, H0, Nnow, z_eq, OmegaT, OmegaL, OmegaK

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
subroutine ion(Ndim, z, Y, f)
    use precision, only : dp
    use constants, only : pi, c, h_P, k_B
    use constants, only : m_1H, not4
    use constants, only : Lambda_H, Lambda_He
    use constants, only : L_He_2p
    use constants, only : A2P_s, A2P_t, sigma_He_2Ps, sigma_He_2Pt, L_He_2Pt, L_He_2St, L_He2St_ion
    use constants, only : CDB_H, CDB_He, CR, CK_H, CK_He, CL_H, CL_He, CT, Bfact
    implicit none

    integer Ndim, Heflag, Heswitch, Hswitch

    real(dp) :: z, x, n, n_He, Trad, Tmat, x_H, x_He
    real(dp) :: y(Ndim), f(Ndim)
    real(dp) :: H_frac
    real(dp) :: Tnow, H0, Nnow, z_eq, Hz, OmegaT, OmegaL, OmegaK
    real(dp) :: Rup, Rdown, K, K_He, Rup_He, Rdown_He, He_Boltz
    real(dp) :: timeTh, timeH, factor
    real(dp) :: fHe, fu, b_He
    real(dp) :: a_VF, b_VF, T_0, T_1, sq_0, sq_1, a_PPB, b_PPB, c_PPB, d_PPB
    real(dp) :: tauHe_s, pHe_s
    real(dp) :: Doppler, gamma_2Ps, pb, qb, AHcon
    real(dp) :: a_trip, b_trip, Rdown_trip, Rup_trip
    real(dp) :: tauHe_t, pHe_t, CL_PSt, CfHe_t, gamma_2Pt
    real(dp) :: AGauss1, AGauss2, zGauss1, zGauss2, wGauss1, wGauss2
    real(dp) :: dHdz, epsilon

    common/Cdata/H_frac, fHe, fu
    common/Hemod/b_He
    common/Hmod/AGauss1, AGauss2, zGauss1, zGauss2, wGauss1, wGauss2
    common/Switch/Heswitch, Hswitch
    common/Cosmo/Tnow, H0, Nnow, z_eq, OmegaT, OmegaL, OmegaK

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
    a_trip = 10._dp**(-16.306_dp)
    b_trip = 0.761_dp

    x_H = y(1)
    x_He = y(2)
    x = x_H + fHe * x_He
    Tmat = y(3)

    n = Nnow * (1._dp + z)**3
    n_He = fHe * Nnow * (1._dp + z)**3
    Trad = Tnow * (1._dp + z)
    Hz = H0 * sqrt(OmegaT * (1._dp + z)**4 / (1._dp + z_eq) &
                   + OmegaT * (1._dp + z)**3 &
                   + OmegaK * (1._dp + z)**2 &
                   + OmegaL)

    ! Also calculate derivative for use later
    dHdz = (H0**2 / 2._dp / Hz) * (4._dp * OmegaT * (1._dp + z)**3 / (1._dp + z_eq) &
                                   + 3._dp * OmegaT * (1._dp + z)**2 &
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
        ! Peebles coefficient K=lambda_a^3/8piH
        K = CK_H / Hz
    else
        ! fit a double Gaussian correction function
        K = CK_H / Hz * (1.0_dp + AGauss1 * exp(-((log(1.0_dp + z) - zGauss1) / wGauss1)**2._dp) &
                                + AGauss2 * exp(-((log(1.0_dp + z) - zGauss2) / wGauss2)**2._dp))
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
    if (Heflag == 0)then        !use Peebles coeff. for He
        K_He = CK_He / Hz
    else    !for Heflag>0       !use Sobolev escape probability
        tauHe_s = A2P_s * CK_He * 3._dp * n_He * (1._dp - x_He) / Hz
        pHe_s = (1._dp - exp(-tauHe_s)) / tauHe_s
        K_He = 1._dp / (A2P_s * pHe_s * 3._dp * n_He * (1._dp - x_He))
        ! smoother criterion here from Antony Lewis & Chad Fendt
        if (((Heflag == 2) .or. (Heflag >= 5)) .and. (x_H < 0.9999999_dp))then
            ! use fitting formula for continuum opacity of H
            ! first get the Doppler width parameter
            Doppler = 2._dp * k_B * Tmat / (m_1H * not4 * c * c)
            Doppler = c * L_He_2p * sqrt(Doppler)
            gamma_2Ps = 3._dp * A2P_s * fHe * (1._dp - x_He) * c * c / (sqrt(pi) * sigma_He_2Ps * 8._dp * pi * Doppler * (1._dp - x_H)) / ((c * L_He_2p)**2._dp)
            pb = 0.36_dp  !value from KIV (2007)
            qb = b_He
            ! calculate AHcon, the value of A * p_(con, H) for H continuum opacity
            AHcon = A2P_s / (1._dp + pb * (gamma_2Ps**qb))
            K_He = 1._dp / ((A2P_s * pHe_s + AHcon) * 3._dp * n_He * (1._dp - x_He))
        end if
        if (Heflag >= 3) then     !include triplet effects
            tauHe_t = A2P_t * n_He * (1._dp - x_He) * 3._dp
            tauHe_t = tauHe_t /(8._dp * pi * Hz * L_He_2Pt**(3._dp))
            pHe_t = (1._dp - exp(-tauHe_t)) / tauHe_t
            CL_PSt = h_P * c * (L_He_2Pt - L_He_2St) / k_B
            if ((Heflag == 3) .or. (Heflag == 5) .or. (x_H > 0.99999_dp)) then
                ! no H cont. effect
                CfHe_t = A2P_t * pHe_t * exp(-CL_PSt / Tmat)
                CfHe_t = CfHe_t / (Rup_trip + CfHe_t)   !"C" factor for triplets
            else                  !include H cont. effect
                Doppler = 2._dp * k_B * Tmat / (m_1H * not4 * c * c)
                Doppler = c * L_He_2Pt * sqrt(Doppler)
                gamma_2Pt = 3._dp * A2P_t * fHe * (1._dp - x_He) * c * c &
                    /(sqrt(pi) * sigma_He_2Pt * 8._dp * pi * Doppler * (1._dp - x_H)) &
                    /((c * L_He_2Pt)**2._dp)
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
    timeH = 2._dp / (3._dp * H0 * (1._dp + z)**1.5)      !Hubble time

    ! calculate the derivatives
    ! turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    ! (clunky, but seems to work)
    if (x_H > 0.99_dp) then         !don't change at all
        f(1) = 0._dp
    !else if ((x_H > 0.98_dp) .and. (Heflag == 0)) then    !don't modify
    else if (x_H > 0.985_dp) then     !use Saha rate for Hydrogen
        f(1) = (x * x_H * n * Rdown - Rup * (1._dp - x_H) * exp(-CL_H / Tmat)) &
            /(Hz * (1._dp + z))
        ! for interest, calculate the correction factor compared to Saha
        ! (without the fudge)
        factor = (1._dp + K * Lambda_H * n * (1._dp - x_H)) &
            /(Hz * (1._dp + z) * (1._dp + K * Lambda_H * n * (1._dp - x) &
            +K * Rup * n * (1._dp - x)))
    else                  !use full rate for H
        f(1) = ((x * x_H * n * Rdown - Rup * (1._dp - x_H) * exp(-CL_H / Tmat)) &
            *(1._dp + K * Lambda_H * n * (1._dp - x_H))) &
            /(Hz * (1._dp + z) * (1._dp / fu + K * Lambda_H * n * (1._dp - x_H) / fu &
            +K * Rup * n * (1._dp - x_H)))
    end if
    ! turn off the He once it is small
    if (x_He < 1.e-15_dp) then
        f(2) = 0._dp
    else
        f(2) = ((x * x_He * n * Rdown_He &
            - Rup_He * (1._dp - x_He) * exp(-CL_He / Tmat)) &
            *(1._dp+ K_He * Lambda_He * n_He * (1._dp - x_He) * He_Boltz)) &
            /(Hz * (1._dp + z) &
            * (1._dp + K_He * (Lambda_He + Rup_He) * n_He * (1._dp - x_He) * He_Boltz))
        ! Modification to HeI recombination including channel via triplets
        if (Heflag >= 3) then
            f(2) = f(2)+ (x * x_He * n * Rdown_trip &
                - (1._dp - x_He) * 3._dp * Rup_trip * exp(-h_P * c * L_He_2St / (k_B * Tmat))) &
                *CfHe_t / (Hz * (1._dp + z))
        end if
    end if

    ! follow the matter temperature once it has a chance of diverging

    if (timeTh < H_frac * timeH) then
        !f(3) = Tmat / (1._dp + z)  !Tmat follows Trad
        ! additional term to smooth transition to Tmat evolution,
        ! (suggested by Adam Moss)
        epsilon = Hz * (1._dp + x + fHe) / (CT * Trad**3 * x)
        f(3) = Tnow &
            + epsilon * ((1._dp + fHe) / (1._dp + fHe + x)) * ((f(1) + fHe * f(2)) / x) &
            - epsilon * dHdz / Hz + 3.0_dp * epsilon / (1._dp + z)
    else
        f(3)= CT * (Trad**4) * x / (1._dp + x + fHe) &
            * (Tmat - Trad) / (Hz * (1._dp + z)) + 2._dp * Tmat / (1._dp + z)
    end if
    return
end subroutine ion


!##############################################################################
      subroutine dverk (n, fcn, x, y, xend, tol, ind, c, nw, w)
        use precision, only : dp
        implicit none

        integer n, ind, nw, k
        real(dp) :: x, y(n), xend, tol, c(24), w(nw,9), temp

        external fcn

!       ******************************************************************
!       * begin initialization, parameter checking, interrupt re-entries *
!       ******************************************************************

!  ......abort if ind out of range 1 to 6
        if (ind < 1 .or. ind > 6) go to 500

        ! cases - initial entry, normal re-entry, interrupt re-entries
        go to (5, 5, 45, 1111, 2222, 2222), ind
        ! case 1 - initial entry (ind == 1 or 2)
!  .........abort if n > nw or tol <= 0
    5   if (n > nw .or. tol <= 0._dp) go to 500
        if (ind == 1) then
            ! initial entry without options (ind == 1)
            ! set c(1) to c(9) equal to 0
            do k = 1, 9
                c(k) = 0._dp
            end do
        else if (ind == 2) then
            ! initial entry with options (ind == 2)
            ! make c(1) to c(9) non-negative
            do k = 1, 9
                c(k) = abs(c(k))
            end do
            ! make floor values non-negative if they are to be used
            if (abs(c(1) - 4._dp) < tiny(1._dp) .or. abs(c(1) - 5._dp) < tiny(1._dp)) then  ! if c(1) == 4 .or. c(1) == 5
                do k = 1, n
                    c(k + 30) = abs(c(k + 30))
                end do
            end if
        end if
        ! initialize rreb, dwarf, prev xend, flag, counts
        c(10) = 2._dp**(-56)
        c(11) = 1.e-35_dp
        ! set previous xend initially to initial value of x
        c(20) = x
        do k = 21, 24
            c(k) = 0._dp
        end do
        go to 50
        ! case 2 - normal re-entry (ind == 3)
!  .........abort if xend reached, and either x changed or xend not
        ! if c(21) /= 0  .and.  ( x /= c(20) .or. xend == c(20) )
   45   if (abs(c(21)) >= tiny(1._dp) .and. (abs(x - c(20)) >= tiny(1._dp) .or. abs(xend - c(20)) < tiny(1._dp))) go to 500
!           re-initialize flag
            c(21) = 0._dp
            ! case 3 - re-entry following an interrupt (ind == 4 to 6)
            ! transfer control to the appropriate re-entry point..........
            ! this has already been handled by the computed go to        .
            ! end cases                                                     v
   50   continue

!       end initialization, etc.

!       ******************************************************************
!       * loop through the following 4 stages, once for each trial  step *
!       * until the occurrence of one of the following                   *
!       *    (a) the normal return (with ind == 3) on reaching xend in *
!       *        stage 4                                                 *
!       *    (b) an error return (with ind < 0) in stage 1 or stage 4 *
!       *    (c) an interrupt return (with ind  ==  4,  5  or  6),  if *
!       *        requested, in stage 1 or stage 4                        *
!       ******************************************************************

99999   continue

!       ***************************************************************
!       * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!       * and some parameter  checking,  and  end  up  with  suitable *
!       * values of hmag, xtrial and htrial in preparation for taking *
!       * an integration step.                                        *
!       ***************************************************************

!***********error return (with ind=-1) if no of fcn evals too great
        if (abs(c(7)) >= tiny(1._dp) .and. c(24) >= c(7)) then  ! if c(7) /= 0 .and. c(24) >= c(7)
            ind = -1
            return
        end if

        ! calculate slope (adding 1 to no of fcn evals) if ind /= 6
        if (ind /= 6) then
            call fcn(n, x, y, w(1,1))
            c(24) = c(24) + 1._dp
        end if

        ! calculate hmin - use default unless value prescribed
        c(13) = c(3)
        if (abs(c(3)) < tiny(1._dp)) then  ! if c(3) == 0
            ! calculate default value of hmin
            ! first calculate weighted norm y - c(12) - as specified
            ! by the error control indicator c(1)
            temp = 0._dp
            if (abs(c(1) - 1._dp) < tiny(1._dp)) then  ! if c(1) == 1
                ! absolute error control - weights are 1
                do k = 1, n
                    temp = max(temp, abs(y(k)))
                end do
                c(12) = temp
            else if (abs(c(1) - 2._dp) < tiny(1._dp)) then  ! if c(1) == 2
                ! relative error control - weights are 1/abs(y(k)) so
                ! weighted norm y is 1
                c(12) = 1._dp
            else if (abs(c(1) - 3._dp) < tiny(1._dp)) then  ! if c(1) == 3
                ! weights are 1/max(c(2), abs(y(k)))
                do k = 1, n
                    temp = max(temp, abs(y(k)) / c(2))
                end do
                c(12) = min(temp, 1._dp)
            else if (abs(c(1) - 4._dp) < tiny(1._dp)) then  ! if c(1) == 4
                ! weights are 1/max(c(k + 30), abs(y(k)))
                do k = 1, n
                    temp = max(temp, abs(y(k)) / c(k + 30))
                end do
                c(12) = min(temp, 1._dp)
            else if (abs(c(1) - 5._dp) < tiny(1._dp)) then  ! if c(1) == 5
                ! weights are 1 / c(k + 30)
                do k = 1, n
                    temp = max(temp, abs(y(k)) / c(k + 30))
                end do
                c(12) = temp
            else
                ! default case - weights are 1/max(1, abs(y(k)))
                do k = 1, n
                    temp = max(temp, abs(y(k)))
                end do
                c(12) = min(temp, 1._dp)
            end if
            c(13) = 10._dp * max(c(11), c(10) * max(c(12) / tol, abs(x)))
        end if

        ! calculate scale - use default unless value prescribed
        c(15) = c(5)
        if (abs(c(5)) < tiny(1._dp)) then  ! if c(5) == 0
            c(15) = 1._dp
        end if

        ! calculate hmax - consider 4 cases
        ! case 1 both hmax and scale prescribed
        if (abs(c(6)) >= tiny(1._dp) .and. abs(c(5)) >= tiny(1._dp)) then  ! if c(6) /= 0 .and. c(5) /= 0
            c(16) = min(c(6), 2._dp / c(5))
        ! case 2 - hmax prescribed, but scale not
        else if (abs(c(6)) >= tiny(1._dp) .and. abs(c(5)) < tiny(1._dp)) then  ! if c(6) /= 0 .and. c(5) == 0
            c(16) = c(6)
        ! case 3 - hmax not prescribed, but scale is
        else if (abs(c(6)) < tiny(1._dp) .and. abs(c(5)) >= tiny(1._dp)) then  ! if c(6) == 0 .and. c(5) /= 0
            c(16) = 2._dp / c(5)
        ! case 4 - neither hmax nor scale is provided
        else if (abs(c(6)) < tiny(1._dp) .and. abs(c(5)) < tiny(1._dp)) then  ! if c(6) == 0 .and. c(5) == 0
            c(16) = 2._dp
        end if

!***********error return (with ind=-2) if hmin > hmax
        if (c(13) > c(16)) then
            ind = -2
            return
        end if

        ! calculate preliminary hmag - consider 3 cases
        if (ind <= 2) then
            ! case 1 - initial entry - use prescribed value of hstart, if
            ! any, else default
            c(14) = c(4)
            if (abs(c(4)) < tiny(1._dp)) c(14) = c(16) * tol**(1. / 6.)  ! if c(4) == 0
        else if (c(23) <= 1._dp) then
            ! case 2 - after a successful step, or at most  one  failure,
            ! use min(2, .9 * (tol/est)**(1/6)) * hmag, but avoid possible
            ! overflow. then avoid reduction by more than half.
            temp = 2._dp * c(14)
                if (tol < (2._dp / .9_dp)**6 * c(19)) &
                    temp = .9_dp * (tol / c(19))**(1. / 6.) * c(14)
            c(14) = max(temp, .5_dp * c(14))
        else
            ! case 3 - after two or more successive failures
            c(14) = .5_dp * c(14)
        end if

        ! check against hmax
        c(14) = min(c(14), c(16))

        ! check against hmin
        c(14) = max(c(14), c(13))

!***********interrupt no 1 (with ind=4) if requested
        if (abs(c(8)) >= tiny(1._dp)) then  ! if c(8) /= 0
            ind = 4
            return
        end if
!       resume here on re-entry with ind == 4   ........re-entry..
 1111   continue

        ! calculate hmag, xtrial - depending on preliminary hmag, xend
        if (c(14) < abs(xend - x)) then
            ! do not step more than half way to xend
            c(14) = min(c(14), .5_dp * abs(xend - x))
            c(17) = x + dsign(c(14), xend - x)
        else
            ! hit xend exactly
            c(14) = abs(xend - x)
            c(17) = xend
        end if

        ! calculate htrial
        c(18) = c(17) - x

!       end stage 1

!       ***************************************************************
!       * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!       * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!       * stage 3. w(*,9) is temporary storage until finally it holds *
!       * ytrial.                                                     *
!       ***************************************************************

        temp = c(18) / 1398169080000._dp

        do k = 1, n
            w(k,9) = y(k) + temp * w(k,1) * 233028180000._dp
        end do
        call fcn(n, x + c(18) / 6._dp, w(1,9), w(1,2))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 74569017600._dp &
                                      + w(k,2) * 298276070400._dp  )
        end do
        call fcn(n, x + c(18) * (4._dp / 15._dp), w(1,9), w(1,3))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 1165140900000._dp &
                - w(k,2) * 3728450880000._dp &
                + w(k,3) * 3495422700000._dp )
        end do
        call fcn(n, x + c(18) * (2._dp / 3._dp), w(1,9), w(1,4))

        do k = 1, n
            w(k,9) = y(k) + temp * ( - w(k,1) * 3604654659375._dp &
                + w(k,2) * 12816549900000._dp &
                - w(k,3) * 9284716546875._dp &
                + w(k,4) * 1237962206250._dp )
        end do
        call fcn(n, x + c(18) * (5._dp / 6._dp), w(1,9), w(1,5))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 3355605792000._dp &
                - w(k,2) * 11185352640000._dp &
                + w(k,3) * 9172628850000._dp &
                - w(k,4) * 427218330000._dp &
                + w(k,5) * 482505408000._dp  )
        end do
        call fcn(n, x + c(18), w(1,9), w(1,6))

        do k = 1, n
            w(k,9) = y(k) + temp * ( - w(k,1) * 770204740536._dp &
                + w(k,2) * 2311639545600._dp &
                - w(k,3) * 1322092233000._dp &
                - w(k,4) * 453006781920._dp &
                + w(k,5) * 326875481856._dp  )
        end do
        call fcn(n, x + c(18) / 15._dp, w(1,9), w(1,7))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 2845924389000._dp &
                - w(k,2) * 9754668000000._dp &
                + w(k,3) * 7897110375000._dp &
                - w(k,4) * 192082660000._dp &
                + w(k,5) * 400298976000._dp &
                + w(k,7) * 201586000000._dp  )
        end do
        call fcn(n, x + c(18), w(1,9), w(1,8))

        ! calculate ytrial, the extrapolated approximation and store
        ! in w(*,9)
        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 104862681000._dp &
                + w(k,3) * 545186250000._dp &
                + w(k,4) * 446637345000._dp &
                + w(k,5) * 188806464000._dp &
                + w(k,7) * 15076875000._dp &
                + w(k,8) * 97599465000._dp   )
        end do

        ! add 7 to the no of fcn evals
        c(24) = c(24) + 7._dp

!       end stage 2

!       ***************************************************************
!       * stage 3 - calculate the error estimate est. first calculate *
!       * the  unweighted  absolute  error  estimate vector (per unit *
!       * step) for the unextrapolated approximation and store it  in *
!       * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!       * specified by the error  control  indicator  c(1).  finally, *
!       * modify  this result to produce est, the error estimate (per *
!       * unit step) for the extrapolated approximation ytrial.       *
!       ***************************************************************

        ! calculate the unweighted absolute error estimate vector
        do k = 1, n
            w(k,2) = (   w(k,1) * 8738556750._dp &
                + w(k,3) * 9735468750._dp &
                - w(k,4) * 9709507500._dp &
                + w(k,5) * 8582112000._dp &
                + w(k,6) * 95329710000._dp &
                - w(k,7) * 15076875000._dp &
                - w(k,8) * 97599465000._dp) / 1398169080000._dp
        end do

        ! calculate the weighted max norm of w(*,2) as specified by
        ! the error control indicator c(1)
        temp = 0._dp
        if (abs(c(1) - 1._dp) < tiny(1._dp)) then  ! if c(1) == 1
            ! absolute error control
            do k = 1, n
                temp = max(temp, abs(w(k,2)))
            end do
        else if (abs(c(1) - 2._dp) < tiny(1._dp)) then  ! if c(1) == 2
            ! relative error control
            do k = 1, n
                temp = max(temp, abs(w(k,2) / y(k)))
            end do
        else if (abs(c(1) - 3._dp) < tiny(1._dp)) then  ! if c(1) == 3
            ! weights are 1/max(c(2), abs(y(k)))
            do k = 1, n
                temp = max(temp, abs(w(k,2)) &
                    / max(c(2), abs(y(k))) )
            end do
        else if (abs(c(1) - 4._dp) < tiny(1._dp)) then  ! if c(1) == 4
            ! weights are 1/max(c(k + 30), abs(y(k)))
            do k = 1, n
                temp = max(temp, abs(w(k,2)) / max(c(k + 30), abs(y(k))) )
            end do
        else if (abs(c(1) - 5._dp) < tiny(1._dp)) then  ! if c(1) == 5
            ! weights are 1/c(k + 30)
            do k = 1, n
                temp = max(temp, abs(w(k,2) / c(k + 30)))
            end do
        else
            ! default case - weights are 1/max(1, abs(y(k)))
            do k = 1, n
                temp = max(temp, abs(w(k,2)) &
                    / max(1._dp, abs(y(k))) )
            end do
        end if

        ! calculate est - (the weighted max norm of w(*,2)) * hmag * scale
        ! - est is intended to be a measure of the error  per  unit
        ! step in ytrial
        c(19) = temp * c(14) * c(15)

!       end stage 3

!       ***************************************************************
!       * stage 4 - make decisions.                                   *
!       ***************************************************************

        ! set ind=5 if step acceptable, else set ind=6
        ind = 5
        if (c(19) > tol) ind = 6

!***********interrupt no 2 if requested
        if (abs(c(9)) >= tiny(1._dp)) then  ! if c(9) /= 0
            return
        end if
!       resume here on re-entry with ind == 5 or 6   ...re-entry..
 2222   continue

        if (ind /= 6) then
            ! step accepted (ind == 5), so update x, y from xtrial,
            ! ytrial, add 1 to the no of successful steps, and set
            ! the no of successive failures to zero
            x = c(17)
            do k = 1, n
                y(k) = w(k,9)
            end do
            c(22) = c(22) + 1._dp
            c(23) = 0._dp
!**************return(with ind=3, xend saved, flag set) if x == xend
            if (abs(x - xend) < tiny(1._dp)) then  ! if x == xend
                ind = 3
                c(20) = xend
                c(21) = 1._dp
                return
            end if
        else
            ! step not accepted (ind == 6), so add 1 to the no of
            ! successive failures
            c(23) = c(23) + 1._dp
!**************error return (with ind=-3) if hmag <= hmin
            if (c(14) <= c(13)) then
                ind = -3
                return
            end if
        end if

!       end stage 4

        go to 99999
!       end loop

!       begin abort action
  500   continue

        write(*,505) ind, tol, x, n, c(13), xend, nw, c(16), c(20), &
            c(22), c(23), c(24), (y(k), k = 1, n)
  505   format( /// 1h0, 58hcomputation stopped in DVERK with the followin &
            g values - &
            / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =, &
                1pd22.15 &
            / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =, &
                1pd22.15 &
            / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =, &
                1pd22.15 &
            / 1h0, 14x, 27hno of successful steps    =, 0pf8.0 &
            / 1h , 14x, 27hno of successive failures =, 0pf8.0 &
            / 1h , 14x, 27hno of function evals      =, 0pf8.0 &
            / 1h0, 23hthe components of y are &
            // (1h , 1p5d24.15)                                           )

!       end abort action

      end subroutine dverk
