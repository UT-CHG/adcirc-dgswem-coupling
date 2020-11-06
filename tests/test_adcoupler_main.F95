module test_ADCouplerMain

    use ADCouplerMain, only: adcoupler_main

    implicit none


contains

    subroutine test_adcoupler_main

        call adcoupler_main()

    end subroutine

end module
