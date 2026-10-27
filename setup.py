from Cython.Build import cythonize
from setuptools import Extension, setup

# No NumPy C-API: arrays are accessed through typed memoryviews (buffer protocol), so the
# extension builds without NumPy headers and works with any NumPy at runtime.
setup(
    ext_modules=cythonize(
        Extension("cpdb.parser", ["cpdb/parser.pyx"]),
        compiler_directives={"language_level": "3"},
    ),
)
