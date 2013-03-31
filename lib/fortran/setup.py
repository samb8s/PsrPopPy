#!/usr/bin/python

"""
setup.py file for SWIG example
"""

from numpy.distutils.core import setup, Extension
import glob

example_module = Extension('_examplecunt',
                          #include_dirs = ['.'],
                           #sources = glob.glob('dmdist.f'),
                           sources=['ykarea.f'],
                           )

setup (name = 'examplecunt',
       version = '0.1',
       author      = "SWIG Docs",
       description = """Simple swig example from docs""",
       ext_modules = [example_module],
       py_modules = ["examplecunt"],
       )
