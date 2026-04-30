from setuptools import setup, Extension
import numpy

module = Extension(
    "waterlib",                  # name of the compiled module
    sources=["waterlib.c"],       # your C file
    include_dirs=[numpy.get_include()],  # include numpy headers
    extra_compile_args=["-O3"],   # optimize!
)

setup(
    name="waterlib",
    ext_modules=[module],
)

