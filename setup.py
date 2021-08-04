from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np


extensions = [
    Extension(
        'pyrecfast',
        ['src/pyrecfast.pyx'],
        extra_link_args = ['build/recfast_wrapper.o']
    )
]

compiler_directives = {'language_level': 3}

setup(
    name='pyrecfast',
    ext_modules=cythonize(extensions, **compiler_directives),
    include_dirs=[np.get_include()],
    define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]
)

