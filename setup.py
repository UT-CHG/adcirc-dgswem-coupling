#!/usr/bin/env python
import os
from setuptools import setup, find_packages

with open(os.path.join(os.path.dirname(__file__),'srcpy','README.md')) as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()


pyadcircdgswem_cmds = ['pyadcircdgswem = pyadcircdgswem.__main__:main']

setup(
    name='pyadcircdgswem',
    version='0.1.0',
    description='The ADCIRC DGSWEM coupler\'s  Python interface for accessing'
                'adcircdgswem\'s variables and subroutines in Python',
    keywords='pyadcircdgswem, Python, interface, f2py',
    long_description=readme,
    author='Gajanan Choudhary',
    author_email='gajananchoudhary91@gmail.com',
    url='https://github.com/gajanan-choudhary/adcirc-dgswem-coupling',
    license=license,
    entry_points={'console_scripts': pyadcircdgswem_cmds},
    package_dir={'': '.'},
    package_data={'': ['libpyadcircdgswem_shared.so', 'pyadcircdgswem*.so']},
    packages=find_packages(exclude=('src', 'srcpy', 'tests', 'doc', 'cmake',
        'build')),
)
