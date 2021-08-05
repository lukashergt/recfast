import sys
import os
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

def test_example_data():
    z_data, x_data = np.loadtxt("test/example_data/example_new_CODATA_AME_2photon_ionisationlevel.out", unpack=True)
    z, x = np.loadtxt("test.out", unpack=True)
    assert_array_equal(z, z_data)
    assert_array_equal(x, x_data)

def test_pyrecfast():
    sys.path.insert(0, "/home/hergtl/Documents/Projects/RecFastPrj/recfast/")
    from pyrecfast import recfast
    
    z_data, x_data = np.loadtxt("test/example_data/example_new_CODATA_AME_2photon_ionisationlevel.out", unpack=True)
    z, x = recfast(OmegaB=0.04, OmegaC=0.20, OmegaL=0.76, H0=70, Tnow=2.725, Yp=0.25, Hswitch=1, Heswitch=6)
    assert_array_equal(z, z_data)
    assert_array_equal(x, x_data)

