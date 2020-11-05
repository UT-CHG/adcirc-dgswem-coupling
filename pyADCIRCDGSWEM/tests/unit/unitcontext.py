#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# pyadcircdgswem - The Python interface of ADCIRC DG-SWEM coupler
# License: BSD 3-Clause License
# Copyright (c) 2020, UT Austin Computational Hydrualics Group
#
"""
Context module for adding system path for unit tests.
"""
import os
import sys
pyadcircdgswem_root=os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__)
            ))))
        )
sys.path.insert(0, os.path.abspath(pyadcircdgswem_root))

#Trigger addition of pyadcircdgswem library path to system
from pyADCIRCDGSWEM import pyadcircdgswem_path

import pyadcircdgswem as pyad
from pyadcircdgswem import pyadmain as pmain
from pyadcircdgswem import utilities as pu

__all__ = ['pyad', 'pmain', 'pu',]

