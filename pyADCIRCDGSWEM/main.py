#!/usr/bin/env python3
#------------------------------------------------------------------------------#
# pyadcircdgswem - The Python interface of ADCIRC DG-SWEM coupler
# License: BSD 3-Clause License
# Copyright (c) 2020, UT Austin Computational Hydrualics Group
#
"""
Main function of pyadcircdgswem: calls initialize, run, and finalize functions.
"""

import time

import pyadcircdgswem_path
import pyadcircdgswem as pyad

################################################################################
DEBUG_LOCAL = 1

__all__ = ['main'] # The only thing from this module to import if needed.

################################################################################
def main():
    """\
    The main function of pyadcircdgswem.

    This function is the equivalent of ADCIRC DG-SWEM coupler's
    adcircdgswem_main.F90 file containing the main program. It calls
    the initialize, run, and finalize functions of the coupler.
    """
    pyad_comm_init = pyad.pyad_comm_module.pyad_comm_init
    pyad_comm_final = pyad.pyad_comm_module.pyad_comm_final
    pyad_read = pyad.pyad_read_module.pyad_read
    pyad_drog_initialize = pyad.pyad_drog.pyad_drog_initialize
    pyad_drog_timestep = pyad.pyad_drog.pyad_drog_timestep

    print("\nRunning ADCIRC DG-SWEM coupler using its Python interface," \
            "pyadcircdgswem\n")

    print("Initializing ADCIRC DG-SWEM coupler")
    t0 = time.time()
    pyad_init()

    t1 = time.time()
    print("Reading ADCIRC DG-SWEM coupler")
    pyad_run()

    t2 = time.time()
    print("Finalizing ADCIRC DG-SWEM coupler")
    pyad_finalize()

    tInit = t1-t0
    tRun  = t2-t1
    tFin  = t3-t2
    tTot  = t3-t0

    print("Initialize time = {0}".format(tInit))
    print("Run time        = {0}".format(tRun))
    print("Finalize time   = {0}".format(tFin))
    print("Total time      = {0}".format(tTot))

    print("\nFinished running pyadcircdgswem")

################################################################################
if __name__ == '__main__':
    main()

