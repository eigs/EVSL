#!/usr/bin/env python
"""
This contains the setup script to build the python interface for EVSL
"""

from distutils.core import setup, Extension

EVSL_MODULE = Extension('_evsl',
                        ['evsl.i'],
                        include_dirs=['../INC'],
                        library_dirs=['../'],
                        libraries=['blas', 'lapack', 'evsl'],
                        extra_link_args=["-lblas", "-llapack", "-levsl"],
                        swig_opts=['-I../INC'])


setup(name='evsl',
      version='0.0.1',
      author="EVSL Team",
      url="http://www-users.cs.umn.edu/~saad/software/EVSL/",
      description="EVSL Python interface""",
      ext_modules=[EVSL_MODULE],
      py_modules=["evsl"],
     )
