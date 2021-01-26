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

#ifdef HAVE_MPI_MOD
    use mpi, only: MPI_COMM_WORLD
#endif
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
#ifndef HAVE_MPI_MOD
        include 'mpif.h'
#endif
        ! Call DGSWEM's init function.
        call dgswem_init()

        ! Call ADCIRC's init function.
        call adcirc_init(MPI_COMM_WORLD)

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
        logical :: no_mpi_finalize_adcirc = .true.

        ! Call ADCIRC's finalize function, but don't actually call mpi_finalize!
        call adcirc_final(no_mpi_finalize_adcirc)

        ! Call DGSWEM's finalize function.
        call dgswem_fin

    end subroutine


    !--------------------------------------------------------------------------!
end module ADCouplerMain

