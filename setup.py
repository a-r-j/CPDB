import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup

#setup(
#    name='cpdb',
#    version='1.0.0',
#    ext_modules=[
#        Extension(
#            'cpdb.parser',
#            sources=['cpdb/parser.pyx'],
#            include_dirs=[numpy.get_include()]
#        ),
#    ],
#)

setup(
    name='cpdb',
    ext_modules=cythonize(["cpdb/parser.pyx"], compiler_directives={'language_level' : "3"}),
    include_dirs=[numpy.get_include()],
)
