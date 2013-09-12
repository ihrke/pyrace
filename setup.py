#
# ----------- setup.py
#
from setuptools import setup, Extension, find_packages
#from setuptools.command.build import build
#from distutils.command.build import build
import numpy as np
import sys, os

MAJOR=0
MINOR=1
MICRO=0
VERSION='%d.%d.%d' % (MAJOR, MINOR, MICRO)


def gsl_get_include():
    if sys.platform == "win32":
        gsl_include = os.getenv('LIB_GSL')
        if gsl_include is None:
            # Environmental variable LIB_GSL not set, use hardcoded path.
            gsl_include = r"c:\Program Files\GnuWin32\include"
        else:
            gsl_include += "/include"
    else:
        gsl_include = os.popen('gsl-config --cflags').read()[2:-1]

    assert gsl_include != '', "Couldn't find gsl. Make sure it's installed and in the path."
    return gsl_include

def gsl_get_library_dir():
    if sys.platform == "win32":
        lib_gsl_dir = os.getenv('LIB_GSL')
        if lib_gsl_dir is None:
            # Environmental variable LIB_GSL not set, use hardcoded path.
            lib_gsl_dir = r"c:\Program Files\GnuWin32\lib"
        else:
            lib_gsl_dir += "/lib"
    else:
        lib_gsl_dir = os.popen('gsl-config --libs').read().split()[0][2:]

    return lib_gsl_dir

def gsl_get_libraries():
    return ['gsl', 'gslcblas']


#
#------------------------------------------------------------------------
#

dist=setup(
    name="pyrace",
    version=VERSION,
    packages=find_packages(),#["pyrace"],
    author="Matthias Mittner",
    author_email="mittner@uva.nl",
    ext_modules=[Extension("_crace",
                           sources=['pyrace/crace.i', 'pyrace/crace.c'],
                           depends=['pyrace/crace.h'],
                           include_dirs=[np.get_include(), gsl_get_include()],
                           libraries=gsl_get_libraries(),
                           library_dirs=[gsl_get_library_dir()],
                           swig_opts=['-modern'])],
    py_modules=['crace'],
    )

