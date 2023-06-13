import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup

setup(
    name='cpdb',
    ext_modules=cythonize(["cpdb/parser.pyx"], compiler_directives={'language_level' : "3"}),
    include_dirs=[numpy.get_include()],
)
