#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# pyadcircdgswem - The Python interface of ADCIRC DG-SWEM coupler
# License: BSD 3-Clause License
# Copyright (c) 2020, UT Austin Computational Hydrualics Group
#
"""
The main pyadcircdgswem module.

Contains:
    (1) The main function imported as "main", and
    (2) The pyadcircdgswem module/library imported as "pyad".

"""

if __name__ == '__main__':
    import pyadcircdgswem_path
    from main import main
else:
    from . import pyadcircdgswem_path
    from .main import main

import pyadcircdgswem as pyad

