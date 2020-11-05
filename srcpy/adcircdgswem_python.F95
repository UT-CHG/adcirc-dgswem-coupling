!------------------------------------------------------------------------------!
! pyadcircdgswem - The Python interface of ADCIRC DG-SWEM coupler
! License: BSD 3-Clause License
! Copyright (c) 2020, UT Austin Computational Hydrualics Group
!
!------------------------------------------------------------------------------!
! File: adcircdgswem_python.F95
! Author: Gajanan K Choudhary, Postdoctoral Fellow
! Location: The University of Texas at Austin
! Date created: 11/04/2020
! Purpose: Generating a Python API of ADCIRC DG-SWEM coupler using 'f2py'
! Dependencies : (1) srcpy/py*.F95 files translating src *.F95 files to f2py
!------------------------------------------------------------------------------!
!
!------------------------------------------------------------------------------!
! Files dependent on the adcircdgswem source code, i.e., src *.F95 files
#include "pyadcircdgswem_main.F95"

!------------------------------------------------------------------------------!
module utilities
    ! The utilities module of Python API of the ADCIRC DG-SWEM coupler
    ! It contains helper variables for the Python side of the adcircdgswem
    ! Python interface. For example, one way to know if DEBUG compiler flag was
    ! passed to the compiler is by looking at the debug variable below on the
    ! Python side.

    implicit none

    logical, parameter :: on = .true.
    logical, parameter :: off = .false.

#ifdef MPI
    logical, parameter :: messg = .true.
#else
    logical, parameter :: messg = .false.
#endif

#ifdef DEBUG
    logical, parameter :: debug = .true.
#else
    logical, parameter :: debug = .false.
#endif


contains

    ! Nothing, for now.

end module utilities

!------------------------------------------------------------------------------!

