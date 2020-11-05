#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# pyadcircdgswem - The Python interface of ADCIRC DG-SWEM coupler
# License: BSD 3-Clause License
# Copyright (c) 2020, UT Austin Computational Hydrualics Group
#
"""
Main function of pyadcircdgswem that gets invoked on the command line.
"""
if __name__ == "__main__":
    import pyadcircdgswem_path
    from main import main
else:
    from . import pyadcircdgswem_path
    from .main import main

################################################################################
if __name__ == '__main__':
    main()

