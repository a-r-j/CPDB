import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup

setup(
    ext_modules=cythonize(
        Extension(
            "cpdb.parser",
            ["cpdb/parser.pyx"],
            include_dirs=[numpy.get_include()],
        ),
        compiler_directives={"language_level": "3"},
    ),
)
