import numpy
from Cython.Build import cythonize
from setuptools import Extension, setup, find_packages

setup(
    name='cpdb-protein',
    version="0.2.0",
    license="MIT",
    author="Arian Jamasb",
    author_email="arian@jamasb.io",
    url="https://github.com/a-r-j/cpdb",
    ext_modules=cythonize(
        ["cpdb/parser.pyx"],
        compiler_directives={'language_level': "3"}
        ),
    include_dirs=[numpy.get_include()],
    setup_requires=["numpy", "cython"],
    install_requires=["numpy", "pandas", "cython"],
    packages=find_packages()
)
