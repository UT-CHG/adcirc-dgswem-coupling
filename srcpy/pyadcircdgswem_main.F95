!------------------------------------------------------------------------------!
! pyadcircdgswem - The Python interface of ADCIRC DG-SWEM coupler
! License: BSD 3-Clause License
! Copyright (c) 2020, UT Austin Computational Hydrualics Group
!
!------------------------------------------------------------------------------!
! File: pyadcircdgswem_main.F95
! Author: Gajanan K Choudhary, Postdoctoral Fellow
! Location: The University of Texas at Austin
! Date created: 11/04/2020
!
!------------------------------------------------------------------------------!
module pyADMain
    use ADCouplerMain

    implicit none


contains


    !--------------------------------------------------------------------------!
    subroutine pyad_main()
        ! Python wrapper of the main function of the ADCIRC DG-SWEM coupler.

        ! Call the main function of the coupler.
        call adcoupler_main()

    end subroutine


    !--------------------------------------------------------------------------!
    subroutine pyad_initialize()
        ! Python wrapper of the initialize function of the ADCIRC DG-SWEM
        ! coupler.

        ! Call the initialize function of the coupler.
        call adcoupler_initialize()

    end subroutine


    !--------------------------------------------------------------------------!
    subroutine pyad_run()
        ! Python wrapper of the run function of the ADCIRC DG-SWEM coupler.

        ! Call the run function of the coupler.
        call adcoupler_run()

    end subroutine


    !--------------------------------------------------------------------------!
    subroutine pyad_finalize()
        ! Python wrapper of the finalize function of the ADCIRC DG-SWEM coupler.
        ! No other instructions should be added after this subroutine is called.

        ! Call the finalize function of the coupler.
        call adcoupler_finalize()

    end subroutine


    !--------------------------------------------------------------------------!
end module pyADMain

