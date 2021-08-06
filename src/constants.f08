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
    real(dp), parameter :: sigma_thomson = 6.6524587321e-29_dp ! Thomson cross section [m^2]

    ! black body radiation constant for u=aT^4:
    real(dp), parameter :: a_rad = 8._dp * pi**5 * k_B**4 / (15._dp * c**3 * h_P**3) ! [J/m^3/K^4]

    ! cosmological abbreviations
    real(dp), parameter :: kappa = 8._dp * pi * G

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
    real(dp), parameter :: CT = (8._dp / 3._dp) * a_rad * sigma_thomson / m_e / c
    real(dp), parameter :: Bfact = h_P * c * (L_He_2p - L_He_2s) / k_B  ! Bfact = exp((E_2p - E_2s) / kT); Extra Boltzmann factor
end module constants
