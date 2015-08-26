#!/usr/bin/env python

from distutils.core import setup

setup(
    name='primer_design',
    version='1.0.0',
    author='Edmund Lau',
    author_email='elau1@students.unimelb.edu.au',
    packages=['source'],
    entry_points={
        'console_scripts': ['primer_design = source.primer_design:main']
    },
    url='https://github.com/edmundlth/primer_design',
    license='LICENSE.txt',
    description=('FIXME'),
    long_description=('FIXME'),
    install_requires=["Biopython"
    ],
)
