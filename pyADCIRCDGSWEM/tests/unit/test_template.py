#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# pyadcircdgswem - The Python interface of ADCIRC DG-SWEM coupler
# License: BSD 3-Clause License
# Copyright (c) 2020, UT Austin Computational Hydrualics Group
#
"""
Unit test template for Python API of ADCIRC DG-SWEM coupler, pyadcircdgswem.
"""
import unittest

if __name__ == '__main__':
    # Executing test directly over command line
    from unitcontext import *
else:
    # Executing test as module
    from .unitcontext import *

################################################################################
LOCALDEBUG = 0

################################################################################
class TestTemplate(unittest.TestCase):
    """Unit test class template."""

    def setUp(self):
        """Initialize the unit test."""
        pass

    def tearDown(self):
        """Finalize the unit test."""
        pass

    def test_main(self):
        """Run the unit test."""
        pass

################################################################################
if __name__ == '__main__':
    unittest.main()

