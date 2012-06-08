#!/usr/bin/env python
# coding: utf-8

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import sys

if sys.version_info < (2, 6):
    print("This module requires Python >= 2.6")
    sys.exit(0)

description="""
A realization of Multifractal Network Generator.

See:
Palla - Vicsek - Lovász: Multifractal network generator, PNAS, 2010
"""

options = dict(
    name = 'mfng',
    version = '0.3',
    description = 'Multifractal network generator',
    long_description = description,
    license = 'BSD License',

    author = 'Arpad Horvath',
    author_email = 'horvath.arpad.szfvar@gmail.com',
    url = 'http://www.arek.uni-obuda.hu/cxnet/doc/html',

    packages = ['mfng'],
    #scripts = ['scripts/igraph'],
    #test_suite = "igraph.test.suite",

    platforms = 'ALL',
    keywords = ['graph', 'network', 'mathematics', 'math', 'graph theory', 'discrete mathematics', 'complex networks'],
    classifiers = [
      'Development Status :: 4 - Beta',
      'Intended Audience :: Education',
      'Intended Audience :: Science/Research',
      'Intended Audience :: Information Technology',
      'Operating System :: OS Independent',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering',
      'Topic :: Scientific/Engineering :: Information Analysis',
      'Topic :: Scientific/Engineering :: Mathematics',
      'Topic :: Scientific/Engineering :: Physics',
      'Topic :: Software Development :: Libraries :: Python Modules',
      'License :: OSI Approved :: BSD License',
    ]
)

setup(**options)
