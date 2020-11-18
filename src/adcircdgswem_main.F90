!------------------------------------------------------------------------------!
! adcircdgswem - The ADCIRC DG-SWEM coupler
! License: BSD 3-Clause License
! Copyright (c) 2020, UT Austin Computational Hydrualics Group
!
!------------------------------------------------------------------------------!
! File: adcircdgswem_main.F95
! Author: Gajanan K Choudhary, Postdoctoral Fellow
! Location: The University of Texas at Austin
! Date created: 11/04/2020
!
!------------------------------------------------------------------------------!
module ADCouplerMain

    implicit none


contains


    !--------------------------------------------------------------------------!
    subroutine adcoupler_main()
        ! Main/top level function of the ADCIRC DG-SWEM coupler.

        ! Initialize the (parallel) coupler.
        call adcoupler_initialize()

        ! Run the (parallel) coupler.
        call adcoupler_run()

        ! Finalize the (parallel) coupler.
        ! No other instructions should be added after the following call.
        call adcoupler_finalize()
    end subroutine


    !--------------------------------------------------------------------------!
    subroutine adcoupler_initialize()
        ! Initialize function of the ADCIRC DG-SWEM coupler.

        ! Do nothing for now.

    end subroutine


    !--------------------------------------------------------------------------!
    subroutine adcoupler_run()
        ! Run function of the ADCIRC DG-SWEM coupler.

        ! Do nothing for now.

    end subroutine


    !--------------------------------------------------------------------------!
    subroutine adcoupler_finalize()
        ! Finalize function of the ADCIRC DG-SWEM coupler.

        ! Do nothing for now.

    end subroutine


    !--------------------------------------------------------------------------!
end module ADCouplerMain

