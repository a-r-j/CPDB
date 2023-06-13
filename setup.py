import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup

setup(
    name='cpdb',
    version="0.0.1",
    license="MIT",
    author="Arian Jamasb",
    author_email="arian@jamasb.io",
    url="https://github.com/a-r-j/cpdb",
    ext_modules=cythonize(["cpdb/parser.pyx"], compiler_directives={'language_level' : "3"}),
    include_dirs=[numpy.get_include()],
    install_requires=["numpy", "pandas"],
)
