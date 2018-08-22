#!/usr/bin/env python

from distutils.core import setup, Extension


dp_module = Extension('_dp',
                           sources=['dp_wrap.c', 'dp.c'],
                           )

setup (name = 'dp',
       version = '0.1',
       author      = "Martin Engler",
       description = "Charge assignment DP",
       ext_modules = [dp_module],
       py_modules = ["dp"],
       )