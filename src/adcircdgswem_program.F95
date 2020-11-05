!------------------------------------------------------------------------------!
! adcircdgswem - The ADCIRC DG-SWEM coupler
! License: BSD 3-Clause License
! Copyright (c) 2020, UT Austin Computational Hydrualics Group
!
!------------------------------------------------------------------------------!
! File: adcircdgswem_program.F95
! Author: Gajanan K Choudhary, Postdoctoral Fellow
! Location: The University of Texas at Austin
! Date created: 11/04/2020
!
!------------------------------------------------------------------------------!
program AdcircDgswemCouplerProgram
    ! The ADCIRC DG-SWEM Coupler Program.

    use ADCouplerMain, only: adcoupler_main

    implicit none

    ! Call the main function.
    call adcoupler_main()

end program

