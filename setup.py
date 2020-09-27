#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import io
import subprocess

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'qsw_mpi'
DESCRIPTION = 'A framework for Quantum Stochastic Walk Simulation.'
URL = 'https://github.com/Edric-Matwiejew/QSW_MPI'
EMAIL = '21469154@student.uwa.edu.au'
AUTHOR = 'Edric Matwiejew'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '1.0.2'

# What packages are required for this module to be exeuted?
REQUIRED = [
        'numpy', 'scipy', 'mpi4py']

# What packages are optional?
EXTRAS = {
        '.qsw support':['h5py'],
        'visulisation':['matplotlib', 'networkx']
        }

EXTENSION = subprocess.run(
        ['python3-config --extension-suffix'],
        shell = True,
        stdout = subprocess.PIPE).stdout.decode('utf-8').strip()

here = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(here, 'README.rst'), encoding = 'utf-8') as f:
    long_description = '\n' + f.read()

setup(
        name = NAME,
        version = VERSION,
        description = DESCRIPTION,
        long_description = long_description,
        long_description_content_type = 'text/markdown',
        author = AUTHOR,
        author_email = EMAIL,
        python_requires = REQUIRES_PYTHON,
        url = URL,
        packages = find_packages(exclude = ['examples/*']),
        package_data = {
            'qsw_mpi': ['*' + EXTENSION],
            },
        install_requires = REQUIRED,
        extras_require = EXTRAS,
        license = 'GPLv3',
        )

