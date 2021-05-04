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
    

    
    use global_dg, only: statim_DG => statim, dt_DG => DTDP, &
    &    rnday_DG => rnday, sz_DG => sz, DG_ntstep => NT, &
    &    nope_DG =>NOPE, neta_DG =>NETA, nvdll_DG => NVDLL, &
    &    nbdv_DG => NBDV, eta_DG=> ETA2, ESBIN1_DG =>ESBIN1, &
    &    ESBIN2_DG =>ESBIN2, ETIME1_DG => ETIME1, ETIME2_DG =>ETIME2, &
    &    ntcysge_DG => NTCYSGE, nspoolge_DG => NSPOOLGE, &
    &    time_a_DG => TIME_A, iths_DG =>ITHS, NSCOUGE_DG => NSCOUGE, &
    &    ntcyfge_DG => NTCYFGE, qnin1_DG =>QNIN1, qnin2_DG =>QNIN2, &
    &    nvel_DG => NVEL, nbou_DG => NBOU, lbcodei_DG => LBCODEI, &
    &    uu2_DG => UU2, vv2_DG => VV2, nbv_DG => NBV, &
    &    QTIME1_DG => QTIME1, QTIME2_DG =>QTIME2, x_DG => X, y_DG => Y, &
    &    nbfr_DG => NBFR, etiminc_DG => ETIMINC, ftiminc_DG => FTIMINC, &
    &    noute_DG => NOUTE, noutv_DG => NOUTV, nffr_DG => NFFR
    use global, only: statim_AD => statim, dt_AD => DT, &
    &    rnday_AD => rnday, sz_AD => sz, AD_ntstep => NT, &
    &    eta_AD => ETA2, ESBIN1_AD =>ESBIN1, &
    &    ESBIN2_AD =>ESBIN2, qnin1_AD =>QNIN1, qnin2_AD =>QNIN2, &
    &    uu2_AD => UU2, vv2_AD => VV2, QTIME2_AD =>QTIME2, &
    &    ftiminc_AD => FTIMINC, nbfr_AD => NBFR, noute_AD => NOUTE, &
    &    noutv_AD => NOUTV
    use boundaries, only: nope_AD =>NOPE, neta_AD =>NETA, &
    &    nvdll_AD => NVDLL, nbdv_AD => NBDV, nvel_AD => NVEL, &
    &    lbcodei_AD => LBCODEI, nbou_AD => NBOU, nbv_AD => NBV, &
    &    nelev_bc_ad => NPEBC
    use gwce, only: ETIME1_AD => ETIME1, ETIME2_AD =>ETIME2, &
    &    etiminc_AD => ETIMINC
    use sizes_dg, only: dg_dir => DIRNAME, dg_indir => INPUTDIR, &
    &    dg_globdir => GLOBALDIR, dg_locdir => LOCALDIR, &
    &    dg_glbindir => GBLINPUTDIR, dg_writeloc => WRITE_LOCAL_FILES, &
    &    dg_rootdir => ROOTDIR
    use dg_read_input_mod
    use sizes, only: ad_locdir => LOCALDIR
    use mesh, only: dp_AD => DP, xcoord_AD => X, ycoord_AD => Y
    use dg, only: num_flow_edges_DG => NFEDS, norm_X_DG => COSNX, &
    &    norm_Y_DG => SINNX, nfedn_DG => NFEDN

#ifdef HAVE_MPI_MOD
    use mpi, only: MPI_COMM_WORLD, MPI_STATUS_SIZE, MPI_INT, &
    &    MPI_REAL, MPI_ANY_SOURCE, MPI_MAX
#endif
    implicit none

    !!!! LOCAL variables!
    real(sz_AD) adc_tprev, adc_tnext, adc_tfinal
    real(sz_DG) dg_tfinal, dg_tprev, dg_tnext
    real(sz_AD) testread, dummyread
    integer i,j,ntsteps,ntsteps_tot_dg
    integer max_neta, max_flow_nvel, num_flux_AD 
    integer num_flux_DG
    integer my_id, ierr, num_procs
    integer status(mpi_status_size)
    integer coupl_DG_id, coupl_AD_id
    integer coupl_DG_id_TEMP, coupl_AD_id_TEMP
    integer global_edge_num
    integer flow_corr_fac
    integer flow_corr_fac_TEMP

    integer unit
    integer etmininc_count
    integer ftmininc_count
    real(sz_AD) time_a_AD

    real,allocatable :: temp_ESBIN1_DG(:),temp_ESBIN2_DG(:)
    real,allocatable :: temp_ESBIN1_AD(:),temp_ESBIN2_AD(:)
    real,allocatable :: temp_QNIN1_DG(:),temp_QNIN2_DG(:)
    real,allocatable :: temp_QNIN1_AD(:),temp_QNIN2_AD(:)


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

        ! We will initially use adcirc to initiate the simulations by starting at time zero.
        ! Note that we set the time step size of DG-SWEM to be half of the ADCIRC time step

        ! We know that the node string in ADCIRC that has specified elevation BCs (thorugh fort.19) is the coupling string.


        call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)
        call mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)



        write(6,*) 'num_proces', num_procs
        write(6,*) 'ierr', ierr
        write(6,*) 'my_id', my_id

        adc_tprev = 0.d0
        dg_tprev = 0.d0
        time_a_AD = 0.d0
        ntsteps = 0
        ntsteps_tot_dg = DG_ntstep
        num_flux_DG = 0
        num_flux_AD = 0
        

        adc_tfinal = (statim_AD+rnday_AD)*86400.00
        dg_tfinal = (statim_DG+rnday_DG)*86400.00

        
        !one adcricstep
        AD_ntstep = 1

        ! Dummy variables used to establish processor id of
        ! the coupling boundaries and sharing of information
        coupl_DG_id_TEMP = -1
        coupl_AD_id_TEMP = -1
        flow_corr_fac_TEMP = -2



        write(6,*) 'elev nodes ADCIRC-', neta_AD
        write(6,*) 'elev nodes DG-SWEM-', neta_DG
        write(6,*) 'nope ADCIRC-', nope_AD
        write(6,*) 'nope DG-SWEM-', nope_DG
        write(6,*) 'nbou ADCIRC-', nbou_AD
        write(6,*) 'nbou DG-SWEM-', nbou_DG
        write(6,*) 'nvel ADCIRC-', nvel_AD
        write(6,*) 'nvel DG-SWEM-', nvel_DG
       

        ! Find elevation BC data
        if(neta_DG.gt.neta_AD) then
            max_neta = neta_DG
        else
            max_neta = neta_AD
        endif


        ! Find flux BC data
        do i = 1, nvel_AD
            if ((lbcodei_AD(i).eq.2).or.(lbcodei_AD(i).eq.12).or.(lbcodei_AD(i).eq.22)) then
                num_flux_AD = num_flux_AD + 1
            endif
        enddo
        !write(6,*) 'num_flux_AD', num_flux_AD
        do i = 1, nvel_DG
            if ((lbcodei_DG(i).eq.2).or.(lbcodei_DG(i).eq.12).or.(lbcodei_DG(i).eq.22)) then
                num_flux_DG = num_flux_DG + 1
            endif
            !write(6,*) 'NBV_DG(i)', nbv_DG(i)
        enddo
        write(6,*) 'num_flux_DG', num_flux_DG

        if(num_flux_DG.gt.num_flux_AD) then
            max_flow_nvel = num_flux_DG
        else
            max_flow_nvel = num_flux_AD
        endif

        !write(6,*) 'max_neta', max_neta
        !write(6,*) 'max_flow_nvel', max_flow_nvel

        
        !allocate these as the max eta
        ALLOCATE (temp_ESBIN1_AD(max_neta),temp_ESBIN2_AD(max_neta))
        ALLOCATE (temp_ESBIN1_DG(max_neta),temp_ESBIN2_DG(max_neta))
        ALLOCATE (temp_QNIN1_AD(max_flow_nvel),temp_QNIN2_AD(max_flow_nvel))
        ALLOCATE (temp_QNIN1_DG(max_flow_nvel),temp_QNIN2_DG(max_flow_nvel))
        
  


        do while (adc_tprev.lt.adc_tfinal)
            
            !Increase step counter
            ntsteps = ntsteps + 1

            !Reset counters and count in time
            etmininc_count = 0
            ftmininc_count = 0
            time_a_AD = dt_AD*ntsteps + statim_AD

            !Since the data structure of the two codes is similar/identical
            !we use this trick to ensure we read from the correct files
            !and correct places in these files
            if ((nbfr_AD.eq.0).and.(nope_AD.gt.0).and.((time_a_AD+dt_AD).gt.ETIME2_AD)) then
                etmininc_count =nint((time_a_AD+dt_AD)/etiminc_AD)*neta_AD + 1
                write(6,*) 'AD etmininc_count:  ', etmininc_count
            endif
            if ((num_flux_AD.gt.0).and.((time_a_AD+dt_AD).gt.QTIME2_AD)) then
                ftmininc_count =nint((time_a_AD+dt_AD)/ftiminc_AD)*neta_AD + 1
                write(6,*) 'AD ftmininc_count:  ', ftmininc_count
            endif

            
            if ((nbfr_AD.eq.0).and.(nope_AD.gt.0).and.(etmininc_count.ne.0)) then
                close(19)
                open(19,file=trim(ad_locdir)//'/'//'fort.19')
                do i = 1 , etmininc_count
                    read(19,*)
                enddo
            endif
            if ((num_flux_AD.gt.0).and.(ftmininc_count.ne.0)) then
                close(20)
                open(20,file=trim(ad_locdir)//'/'//'fort.20')
                do i = 1 , ftmininc_count
                    read(20,*)
                enddo
            endif

 
            
            !Find the processors that contain the interface in each code
            if((ntsteps.eq.1).and.(nbfr_AD.eq.0).and.(nope_AD.gt.0).and.(ESBIN1_AD(1).eq.0.d0)) then
                coupl_AD_id_TEMP = my_id
                !write(6,*) 'coupl_AD_id inside if ', coupl_AD_id_TEMP !! coupl_AD_id = 0
            endif
            !call  mpi_bcast(coupl_AD_id, 1, mpi_int, coupl_AD_id, MPI_COMM_WORLD, ierr)
            if((ntsteps.eq.1).and.(num_flux_DG.gt.0).and.(QNIN1_DG(1).eq.0.d0)) then
                coupl_DG_id_TEMP = my_id
                !write(6,*) 'coupl_DG_id inside if ', coupl_DG_id_TEMP
            endif
            if(ntsteps.eq.1) then
                call mpi_allreduce(coupl_DG_id_TEMP, coupl_DG_id, 1, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
                call mpi_allreduce(coupl_AD_id_TEMP, coupl_AD_id, 1, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
            endif


            !write(6,*) 'my_id,  coupl_DG_id', my_id, coupl_DG_id
            !write(6,*) 'my_id,  coupl_AD_id', my_id, coupl_AD_id

            ! Find normal vector on DG side to determine normal vector
            ! in the ADCIRC domain and use it to ensure
            ! correct flow direction over the coupling boundary
            if((ntsteps.eq.1).and.(my_id.eq.coupl_DG_id)) then
                global_edge_num = nfedn_DG(2)
                !write(6,*) 'global flow edge num DG ', global_edge_num
                !write(6,*) 'NX ', norm_X_DG(global_edge_num)
                !write(6,*) 'NY', norm_Y_DG(global_edge_num)
                if(norm_X_DG(global_edge_num).gt.0) then
                    flow_corr_fac_TEMP = -1
                endif
                if(norm_X_DG(global_edge_num).lt.0) then
                    flow_corr_fac_TEMP = 1
                endif
            endif
            if(ntsteps.eq.1) then
                call mpi_allreduce(flow_corr_fac_TEMP, flow_corr_fac, 1, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
            endif
            !write(6,*) 'flow_corr_fac ', flow_corr_fac


            !Call ADCIRC to do one step in time
            call adcirc_run(AD_ntstep)
            adc_tprev = adc_tprev + dt_AD

            if((nope_AD.gt.0).and.(my_id.eq.coupl_AD_id)) then

                do i = 1, neta_AD
                    !multiply veolcity by depth to get correct
                    !inflow BC, the correction factor ensures
                    temp_QNIN1_DG(i) = uu2_AD(nbdv_AD(nope_AD,i))*flow_corr_fac*(temp_ESBIN2_DG(i)+dp_AD(nbdv_AD(nope_AD,i)))
                    temp_QNIN2_DG(i) = uu2_AD(nbdv_AD(nope_AD,i))*flow_corr_fac*(temp_ESBIN2_DG(i)+dp_AD(nbdv_AD(nope_AD,i)))
                    !write(6,*) 'temp_QNIN2_DG', temp_QNIN2_DG(i)
                enddo
                !Now need to put ADCIRC data into DG code.
                !In case we are on the same processor set the
                !variables directly
                if(nope_AD.eq.num_flux_DG) then
                    j = 1
                    do i = 1, nvel_DG
                        if ((lbcodei_DG(i).eq.2).or.(lbcodei_DG(i).eq.12).or.(lbcodei_DG(i).eq.22)) then
                            QNIN2_DG(i) = temp_QNIN2_DG(j)
                            QNIN1_DG(i) = temp_QNIN2_DG(j)
                            j = j + 1
                        endif
                    enddo
                    QTIME2_DG = QTIME2_DG + dt_DG
                endif
                !If we are not on the same processor we send/recieve
                if((my_id.eq.coupl_AD_id)) then
                    call mpi_send( temp_QNIN2_DG, neta_AD ,mpi_real, coupl_DG_id,10,MPI_COMM_WORLD,ierr)
                endif
            endif
            
            if((num_flux_DG.gt.0).and.(my_id.eq.coupl_DG_id)) then
                call mpi_recv(temp_QNIN2_DG, max_flow_nvel, mpi_real, coupl_AD_id, 10, MPI_COMM_WORLD, status, ierr)
                j = 1
                do i = 1, nvel_DG
                    if ((lbcodei_DG(i).eq.2).or.(lbcodei_DG(i).eq.12).or.(lbcodei_DG(i).eq.22)) then
                        QNIN2_DG(i) = temp_QNIN2_DG(j)
                        QNIN1_DG(i) = temp_QNIN2_DG(j)
                        j = j + 1
                        !write(6,*) 'rec QNIN2_DG(j) ', temp_QNIN2_DG(j-1)
                    endif
                enddo
                QTIME2_DG = QTIME2_DG + dt_DG
            endif

            !Need to find correct spot in .19 and .20 files
            !if these were openend and closed
            etmininc_count = 0
            ftmininc_count = 0
            if ((nbfr_DG.eq.0).and.(nope_DG.gt.0).and.((time_a_DG+dt_DG).gt.ETIME2_DG)) then
                etmininc_count =nint((time_a_DG+dt_DG)/etiminc_DG)*neta_DG + 1
                !write(6,*) 'etmininc_count:  ', etmininc_count
            endif
            if ((num_flux_DG.gt.0).and.((time_a_DG+dt_DG).gt.QTIME2_DG)) then !!!here!!!
                ftmininc_count =nint((time_a_DG+dt_DG)/ftiminc_DG)*neta_DG + 1
                !write(6,*) 'ftmininc_count:  ', ftmininc_count
            endif

            
            ! need to make sure we have the correct file open if we have elevation bc
            if ((nbfr_DG.eq.0).and.(nope_DG.gt.0).and.(etmininc_count.ne.0)) then
                close(19)
                open(19,file=dg_dir//'/'//'fort.19')
                do i = 1 , etmininc_count
                    read(19,*)
                enddo
            endif
            if ((num_flux_DG.gt.0).and.(ftmininc_count.ne.0)) then
                close(20)
                open(20,file=dg_dir//'/'//'fort.20')
                do i = 1 , ftmininc_count
                    read(20,*)
                enddo
            endif
            !At the first time step we need to make these zero
            if (time_a_DG.eq.0) then
                DG_ntstep = 0
                iths_DG = 0
            endif

            dg_tfinal = dg_tprev + dt_DG
            DG_ntstep = DG_ntstep+1

            !Need to check if DG-SWEM has closed files unexpectedly
            inquire(file=dg_dir//'/'//'fort.63', number=unit)
            !If this is the case, reopen to append to these files
            if (unit.eq.-1) then
                if (noute_DG.gt.0) then
                open(61,file=dg_dir//'/'//'fort.61', position="append")
                endif
                if (noutv_DG.gt.0) then
                open(62,file=dg_dir//'/'//'fort.62', position="append")
                endif
                open(63,file=dg_dir//'/'//'fort.63', position="append")
                open(64,file=dg_dir//'/'//'fort.64', position="append")
            endif
            !inquire(file=dg_dir//'/'//'fort.19', number=unit)
            !write(6,*) 'unit 19:  ', unit


            !Output for debugging
            !write(6,*) 'dg_dir:  ', dg_dir
            !write(6,*) 'dg_rootdir: ', trim(dg_rootdir)
            !write(6,*) 'INPUTDIR:  ', trim(dg_indir)
            !write(6,*) 'GLOBALDIR:  ', trim(dg_globdir)
            !write(6,*) 'LOCALDIR:  ', trim(dg_locdir)
            !write(6,*) 'GBLINPUTDIR:  ', trim(dg_glbindir)
            !write(6,*) 'AD LOCALDIR:  ', trim(ad_locdir)

            call dgswem_run



            !trick DG-SWEM to only step once at each call
            dg_tprev = dg_tprev + dt_DG
            iths_DG = iths_DG + 1

            if((num_flux_DG.gt.0).and.(my_id.eq.coupl_DG_id)) then
                j = 1
                do i = 1, nvel_DG
                    if ((lbcodei_DG(i).eq.2).or.(lbcodei_DG(i).eq.12).or.(lbcodei_DG(i).eq.22)) then
                        temp_ESBIN2_AD(j) = eta_DG(nbv_DG(i))
                        temp_ESBIN1_AD(j) = eta_DG(nbv_DG(i))
                        j = j + 1
                        !write(6,*) 'temp_ESBIN2_AD(j)', temp_ESBIN2_AD(j-1)
                    endif
                enddo
                !if we are on the same processor
                if(neta_AD.eq.num_flux_DG) then
                    do i = 1, neta_DG
                        ESBIN2_AD(i) = temp_ESBIN2_AD(i)
                        ESBIN1_AD(i) = temp_ESBIN1_AD(i)
                    enddo
                    ETIME2_AD = ETIME2_AD + dt_AD
                endif
                if((my_id.eq.coupl_DG_id)) then
                    call mpi_send( temp_ESBIN2_AD, num_flux_DG ,mpi_real, coupl_AD_id,10,MPI_COMM_WORLD,ierr)
                endif
            endif
            !if we are on separate processors
            if((nope_AD.gt.0).and.(my_id.eq.coupl_AD_id)) then
                call mpi_recv(temp_ESBIN2_AD,neta_AD, mpi_real, coupl_DG_id, 10, MPI_COMM_WORLD, status, ierr)
                do i = 1, neta_AD
                    ESBIN2_AD(i) = temp_ESBIN2_AD(i)
                    ESBIN1_AD(i) = temp_ESBIN2_AD(i)
                    !write(6,*) 'rec temp_ESBIN2_AD', temp_ESBIN2_AD(i))
                enddo
                ETIME2_AD = ETIME2_AD + dt_AD
            endif

        enddo


        write(6,*) 'dt_AD', dt_AD
        write(6,*) 'dt_DG', dt_DG
        !write(6,*) 'AD_ntstep', AD_ntstep
        !write(6,*) 'DG_ntstep', DG_ntstep

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

