# Python wrapper for RECFAST (recfast.for)

## Dependencies

* setuptools
* Cython
* numpy

## Installation

Ideally it should be as simple as:
```bash
git clone https://github.com/lukashergt/recfast.git
cd recfast
make
```

For the moment this is probably bound to fail for any system other than Linux.

You can test whether things are working by running:
```bash
python -c "from pyrecfast import recfast; print(recfast(0.04, 0.20, 0.76, 70, 2.725, 0.25))"
```



# Archive

The following is taken directly from the original website of Prof. Douglas
Scott (https://www.astro.ubc.ca/people/scott/recfast.html).

## RECFAST

### A code to calculate the recombination history of the Universe

This code is meant to reproduce the calculations described in the paper "A new calculation of the recombination epoch", Seager, S., Sasselov, D. & Scott, D., 1999, ApJ, 523, L1 ([astro-ph/9909275](http://arxiv.org/abs/astro-ph/9909275)), which is a fast approximation to the detailed calculations described in "How exactly did the Universe become neutral?" Seager, S., Sasselov, D. & Scott, D., 2000, ApJS, 128, 407 ([astro-ph/9912182](http://arxiv.org/abs/astro-ph/9912182)). And with some updates (to He etc.) which are described in Wong, W.Y., Moss, A. & Scott, D., 2008, MNRAS, 386, 1023 ([arXiv:0711.1357](http://arxiv.org/abs/0711.1357)), together with the Compton coupling treatment in Scott, D. & Moss A., 2009, MNRAS, 397, 445 ([arXiv:0902.3438](http://arxiv.org/abs/0902.3438)).

It solves a modified 3-level atom for each of hydrogen and helium, along with some corrections motivated by more detailed calculations. It is thus *fast*, while being accurate enough for current CMB experiments.

Here are things that the code treats:

* Solves H and He simultaneously
* Kinetic temperature of matter evolved separately, including Compton and adiabatic cooling
* Careful consideration of dependency of quantities on $T_M$ and $T_R$
* Accurate values of all relevant physical constants
* Full treatment of background cosmology, including radiation, Lambda and curvature (but no explicit "$w$")
* Saha assumed for ionized helium
* Accurate look-up table for recombination coefficients
* Approximation for H recombination at low redshift, to account for out of equilibrium
* Singlet and triplet states considered for neutral He
* Additional "fudge factor" for He recombination
* Solves for temperature difference between matter and radiation, giving smoother transition in $T_M$ evolution
* Includes an extra "fudge" function to correct $K(z)$ to approximate the H physics corrections from detailed codes (see Rubino-Martin et al. 2010)
* Extra "fudge factor" Gaussians, designed to make correction for additional He physics for cosmologies near the standard one

This results in $x_e(z)$ producing (it is believed) better than 0.1 percent accuracy for CMB anisotropy $C_\ell$.


### Original Readme

These are instructions for the code recfast.for, which calculates the ionized
fraction as a function of redshift, as described in Seager, Sasselov & Scott
(1999, ApJ, 523, L1), and in more detail in Seager, Sasselov & Scott (2000,
ApJS, 128, 407), with modifications described in Wong, Moss & Scott (2008).
An additional modification to the matter temperature evolution comes from
working out the approximate difference between Trad and Tmat, and then
including the additional derivative terms for dTmat/dz (essentially described
in Peebles (1971) book and in Appendix A of Hirata (2008) - see Scott &
Moss 2009).  Following calculation of new physical effects on hydrogen
in papers by Chluba et al. and Hirata et al. (collated in Rubino-Martin et
al. 2010), a "fudge function" has been fitted in redshift within recfast,
to give approximately the same x_e.  This specifically a change to the old H
fudge factor, together with a double Gaussian function.  The result of this
is that recfast should produce C_l which are accurate to about the 0.1%
level, for all physical effects known at the beginning of 2010.  A minor
modification of the parameters of the second Gaussian gives better
agreement with the additional He corrections described in Chluba, Fung &
Switzer (2012).  RECFAST should now give results which agree with
COSMOREC and HYREC at about the same level with which they agree with each
other - at least for cosmologies close to the standard one.

The code is written in (fairly) standard FORTRAN77.  [One non-standard
feature is the use of tab spacing at the beginnings of lines, which will
typically not be copied using cut and paste, and tabs for line continuation,
which will cause problems for some compilers.]

It will generally be more efficient to compile the code using "- fast" or
the equivalent.

The comments are fairly extensive and the argument names self-explanatory.
Most of the arguments are described in detail in the comments of the code.
Suggestions for additional comments should be sent to the author.

For simplicity the code is presented as a single unit.  It contains a main
program, the subroutine GET_INIT (to initialise the ionized fraction for
arbitrary starting redshift), the subroutine ION (which calculates the
derivatives) and the numerical ODE solver DVERK.

The input interface was designed to look familiar to users of Seljak &
Zaldarriaga's code CMBFAST.  A convenient way to run the program is by using
a file recfast.run of the form:

```
output.file
Omega_B, Omega_DM, Omega_vac
H_0, T_0, Y_p
Hswitch
Heswitch
```

For example:
```
junk.out
0.04 0.20 0.76
70 2.725 0.25
1
6
```

The generic output is a 2 column fomatted file containing redshift (not 1+z)
and `x_e` (defined as `n_e/n_H`).  The default output is steps of 10 in `z` from
10,000 to 0 (but this is easily changed).

Users wishing to include the code as a subroutine (e.g. in a CMB Boltzmann
code) should find it relatively straightforward to make the required changes.
Similarly it should be easy to add simple variations (to the background
cosmology, for example).

Please cite Seager, Sasselov & Scott (1999, ApJ, 523, L1; or 2000,
ApJS, 128, 407) and/or Wong, Moss & Scott (2008, MNRAS, 386, 1023)
when using this code as part of any published work.
