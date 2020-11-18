!------------------------------------------------------------------------------!
! adcircdgswem - The ADCIRC DG-SWEM coupler
! License: BSD 3-Clause License
! Copyright (c) 2020, UT Austin Computational Hydrualics Group
!
!------------------------------------------------------------------------------!
! File: adcircdgswem_main.F90
! Author: Gajanan K Choudhary, Postdoctoral Fellow
! Location: The University of Texas at Austin
! Date created: 11/04/2020
!
!------------------------------------------------------------------------------!
module ADCouplerMain
    use adcirc_export, only: adcirc_init, adcirc_run, adcirc_final
    use dgswem_mod, only: dgswem_init, dgswem_run, dgswem_fin

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

        ! Call ADCIRC's init function.
        call adcirc_init

        ! Call DGSWEM's init function.
        call dgswem_init

    end subroutine


    !--------------------------------------------------------------------------!
    subroutine adcoupler_run()
        ! Run function of the ADCIRC DG-SWEM coupler.

        ! Do nothing for now.
        ! Time loop that calls ADCIRC and DG-SWEM's run functions and
        ! their modifies boundary conditions.

        ! Just testing for now
        call adcirc_run
        call dgswem_run

    end subroutine


    !--------------------------------------------------------------------------!
    subroutine adcoupler_finalize()
        ! Finalize function of the ADCIRC DG-SWEM coupler.

        ! Call DGSWEM's finalize function.
        call dgswem_fin

        ! Call ADCIRC's finalize function.
        call adcirc_final

    end subroutine


    !--------------------------------------------------------------------------!
end module ADCouplerMain

