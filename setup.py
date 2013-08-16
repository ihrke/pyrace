#
# ----------- setup.py
#
from setuptools import setup, Extension
import numpy as np

setup(
    name="pyrace",
    packages=["pyrace"],
    author="Matthias Mittner",
    author_email="mittner@uva.nl",
    ext_modules=[Extension("_pyrace",
                           sources=['pyrace.i, pyrace.c'],
                           depends=['pyrace.h'],
                           include_dirs=[np.get_include()],
                           swig_opts=['-modern'])]
    )


