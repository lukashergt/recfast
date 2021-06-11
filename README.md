# Python wrapper for RECFAST (recfast.for)

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

output.file
Omega_B, Omega_DM, Omega_vac
H_0, T_0, Y_p
Hswitch
Heswitch

For example:
junk.out
0.04 0.20 0.76
70 2.725 0.25
1
6

The generic output is a 2 column fomatted file containing redshift (not 1+z)
and x_e (defined as n_e/n_H).  The default output is steps of 10 in z from
10,000 to 0 (but this is easily changed).

Users wishing to include the code as a subroutine (e.g. in a CMB Boltzmann
code) should find it relatively straightforward to make the required changes.
Similarly it should be easy to add simple variations (to the background
cosmology, for example).

Please cite Seager, Sasselov & Scott (1999, ApJ, 523, L1; or 2000,
ApJS, 128, 407) and/or Wong, Moss & Scott (2008, MNRAS, 386, 1023)
when using this code as part of any published work.
