#!/usr/bin/env python

from distutils.core import setup

setup(
    name='primer_design',
    version='0.1.0',
    author='Edmund Lau',
    author_email='elau1@student.unimelb.edu.au',
    packages=['primer_design'],
    entry_points={
        'console_scripts': ['primer_design = primer_design.primer_desig3:main']
    },
    url='https://github.com/edmundlth/primer_design',
    licence = "LICENCE.txt",
    description = ('README'),
    long_description = ('README'),
)
