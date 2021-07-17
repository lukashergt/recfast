import numpy as np
from numpy cimport ndarray

cdef extern:
    void c_recfast(double *OmegaB,
                   double *OmegaC,
                   double *OmegaL,
                   double *H0inp,
                   double *Tnow,
                   double *Yp,
                   int *Hswitch_in,
                   int *Heswitch_in,
                   double *zinitial,
                   double *zfinal,
                   double *tol,
                   int *Nz,
                   double *z_array,
                   double *x_array)

def recfast(double OmegaB,
            double OmegaC,
            double OmegaL,
            double H0inp,
            double Tnow,
            double Yp,
            int Hswitch_in,
            int Heswitch_in,
            double zinitial,
            double zfinal,
            double tol,
            int Nz):
    cdef:
        ndarray[double, mode="c"] z_array = np.empty(Nz, dtype=np.double)
        ndarray[double, mode="c"] x_array = np.empty(Nz, dtype=np.double)
    c_recfast(&OmegaB, &OmegaC, &OmegaL, &H0inp, &Tnow, &Yp, &Hswitch_in, &Heswitch_in,
              &zinitial, &zfinal, &tol, &Nz, &z_array[0], &x_array[0])
    return z_array, x_array
