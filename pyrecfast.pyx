import numpy as np
from numpy cimport ndarray

cdef extern:
    void c_recfast(double *OmegaB,
                   double *OmegaC,
                   double *OmegaL,
                   double *H0,
                   double *Tnow,
                   double *Yp,
                   int *Hswitch,
                   int *Heswitch,
                   double *zinitial,
                   double *zfinal,
                   double *tol,
                   int *Nz,
                   double *z_array,
                   double *x_array)

def recfast(double OmegaB,
            double OmegaC,
            double OmegaL,
            double H0,
            double Tnow,
            double Yp,
            int Hswitch=1,
            int Heswitch=6,
            double zinitial=10000,
            double zfinal=0,
            double tol=1e-5,
            int Nz=1000):
    """
    Recfast: Integrator for Cosmic Recombination of Hydrogen and Helium.

    Calculate ionised fraction as a function of redshift. Solves for H and He
    simultaneously, and includes H "fudge factor" for low z effect, as well as
    HeI fudge factor.


    Parameters
    ----------
        OmegaB: float
            Baryon density parameter.

        OmegaC: float
            Cold dark matter density parameter.

        OmegaL: float
            Dark energy density parameter.

        H0: float
            Hubble parameter in km/s/Mpc.

        Tnow: float
            Present-day CMB temperature in K.

        Yp: float
            Present-day Helium fraction.

        Hswitch: int
            Integer switch for modifying the H recombination.
                (0) no change from old Recfast, fudge factor `fu=1.14`
                (1) include correction,         fudge factor `fu=1.125`
            default: 1

        Heswitch: int
            Integer switch for modifying the HeI recombination.
                (0) use Peebles coeff. for He
                    (no change from old Recfast)
                (1) for Heswitch > (0) use Sobolev escape probability for singlet 1P-1S transition
                (2) include Heswitch (1) and use fitting formula for continuum opacity of H on HeI
                    (based on fitting formula by Kholupenko, Ivanchik & Varshalovich, 2007)
                (3) include Heswitch (1) and include triplet effects
                (4) include Heswitch (1) and (3) and include H continuum effects
                (5) include Heswitch (1) to (3)
                (6: include Heswitch (1) to (4)
            default: 6

        zinitial: float
            Initial redshift value.
            default: 10000

        zfinal: float
            Final redshift value.
            default: 0 (today)

        tol: float
            Tolerance for numerical integrator.
            default: 1e-5

        Nz: int
            Number or redshift bins between `zinitial` and `zfinal`.
            Determines the size of the returned arrays.
            default: 1000


    Returns
    -------
        z_array: np.ndarray
            Array of redshift values `z`.
            Same shape as `x_array`.

        x_array: np.ndarray
            Array of ionisation fraction `x_e(z)`.
            Same shape as `z_array`.
    """
    cdef:
        ndarray[double, mode="c"] z_array = np.empty(Nz, dtype=np.double)
        ndarray[double, mode="c"] x_array = np.empty(Nz, dtype=np.double)
    c_recfast(&OmegaB, &OmegaC, &OmegaL, &H0, &Tnow, &Yp, &Hswitch, &Heswitch,
              &zinitial, &zfinal, &tol, &Nz, &z_array[0], &x_array[0])
    return z_array, x_array
