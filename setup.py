from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np


extensions = [
    Extension(
        'pyrecfast',
        ['pyrecfast.pyx'],
        extra_link_args = ['recfast_wrapper.o']
    )
]

compiler_directives = {'language_level': 3}

setup(
    name='pyrecfast',
    ext_modules=cythonize(extensions, **compiler_directives),
    include_dirs=[np.get_include()],
    define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]
)

# setup(
#     name='pyrecfast',
#     ext_modules=cythonize([Extension("pyrecfast", ["pyrecfast.pyx"], extra_link_args=["pyrecfast.o"])],
#                           compiler_directives={'language_level': "3"}),
#     zip_safe=False,
#     include_dirs=[np.get_include()]
# )

