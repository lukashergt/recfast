import numpy as np
from numpy.testing import assert_array_equal


data = "test/example_data/example_new_CODATA_AME_2photon_ionisationlevel.out"


def test_recfast_out_to_example_data():
    z_data, x_data = np.loadtxt(data, unpack=True)
    z, x = np.loadtxt("test/example_data/test.out", unpack=True)
    assert_array_equal(z, z_data)
    assert_array_equal(x, x_data)


def test_pyrecfast_to_example_data():
    from pyrecfast import recfast
    z_data, x_data = np.loadtxt(data, unpack=True)
    z, x = recfast(OmegaB=0.04,
                   OmegaC=0.20,
                   OmegaL=0.76,
                   H0=70,
                   Tnow=2.725,
                   Yp=0.25,
                   Hswitch=1,
                   Heswitch=6)
    assert z.size == x.size == 1000
    assert z[-1] == 0
    assert x[0] > 1
    assert x[-1] < 1e-3
    assert_array_equal(z, z_data)
    assert_array_equal(x, x_data)
