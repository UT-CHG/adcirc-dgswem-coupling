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
    &    noute_DG => NOUTE, noutv_DG => NOUTV, nffr_DG => NFFR, &
    &    ntcyfe_DG => NTCYFE, ntcyfv_DG => NTCYFV, ntcyfgv_DG => NTCYFGV, &
    &    toutfge_DG => TOUTFGE
    use global, only: statim_AD => statim, dt_AD => DT, &
    &    rnday_AD => rnday, sz_AD => sz, AD_ntstep => NT, &
    &    eta_AD => ETA2, ESBIN1_AD =>ESBIN1, &
    &    ESBIN2_AD =>ESBIN2, qnin1_AD =>QNIN1, qnin2_AD =>QNIN2, &
    &    uu2_AD => UU2, vv2_AD => VV2, QTIME2_AD =>QTIME2, &
    &    ftiminc_AD => FTIMINC, nbfr_AD => NBFR, noute_AD => NOUTE, &
    &    noutv_AD => NOUTV, emo_AD => EMO
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
    integer num_neta19_AD
    integer my_id, ierr, num_procs, proc_id
    integer status(mpi_status_size)
    integer global_edge_num
    integer flow_corr_fac
    integer flow_corr_fac_TEMP

    integer unit
    integer etmininc_count
    integer ftmininc_count
    real(sz_AD) time_a_AD
    real(sz_AD) temp_NX, temp_NY
    real(sz_AD) avg_NX, avg_NY
    integer coupl_dir, coupl_dir_TEMP
    integer procsum
    integer nope_AD_coupl_ID

    real,allocatable :: temp_ESBIN1_DG(:),temp_ESBIN2_DG(:)
    real,allocatable :: temp_ESBIN1_AD(:),temp_ESBIN2_AD(:)
    real,allocatable :: temp_QNIN1_DG(:),temp_QNIN2_DG(:)
    real,allocatable :: temp_QNIN1_AD(:),temp_QNIN2_AD(:)
    integer,allocatable ::  coupl_DG_id(:), coupl_AD_id(:)
    integer,allocatable ::  coupl_DG_id_TEMP(:), coupl_AD_id_TEMP(:)
    integer,allocatable ::  coupl_DG_fluxnode(:), coupl_AD_Netanode(:)
    integer,allocatable ::  coupl_DG_fluxnode_T(:), coupl_AD_Netanode_T(:)

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
        ! Time loop that calls ADCIRC and DG-SWEM's run functions and
        ! their modifies boundary conditions.

        ! We will initially use adcirc to initiate the simulations by starting at time zero.
        ! Note that we set the time step size of DG-SWEM to be half of the ADCIRC time step

        ! We know that the node string in ADCIRC that has specified elevation BCs (thorugh fort.19) is the coupling string.


        call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)
        call mpi_comm_size(MPI_COMM_WORLD,num_procs,ierr)

        !Add 1 to my_id to conform
        !my_id = my_id + 1
        write(6,*) 'num_procs', num_procs
        write(6,*) 'ierr', ierr
        write(6,*) 'my_id', my_id

        adc_tprev = 0.d0
        dg_tprev = 0.d0
        time_a_AD = 0.d0
        ntsteps = 0
        ntsteps_tot_dg = DG_ntstep
        num_flux_DG = 0
        num_flux_AD = 0
        procsum = 0
        num_neta19_AD = 0
        

        adc_tfinal = (statim_AD+rnday_AD)*86400.00
        dg_tfinal = (statim_DG+rnday_DG)*86400.00

        
        !one adcricstep
        AD_ntstep = 1


        write(6,*) 'elev nodes ADCIRC-', neta_AD
        write(6,*) 'elev nodes DG-SWEM-', neta_DG
        write(6,*) 'nope ADCIRC-', nope_AD
        write(6,*) 'nope DG-SWEM-', nope_DG
        write(6,*) 'nbou ADCIRC-', nbou_AD
        write(6,*) 'nbou DG-SWEM-', nbou_DG
        write(6,*) 'nvel ADCIRC-', nvel_AD
        write(6,*) 'nvel DG-SWEM-', nvel_DG
       



        if (neta_AD.gt.0) then
            do i = 1, neta_AD
                if ((emo_AD(1,i).eq.0).and.(nbfr_AD.gt.0)) then
                    num_neta19_AD = num_neta19_AD +1
                endif
                if (nbfr_AD.eq.0) then
                    num_neta19_AD = neta_AD
                endif
            enddo
            write(6,*) 'my_id  num_neta19_AD ', my_id, num_neta19_AD
        endif
        if (nope_AD.gt.1) then
            do i = 1, nope_AD
                if (nvdll_AD(i).eq.num_neta19_AD) then
                    nope_AD_coupl_ID = i
                    !write(6,*) 'my_id nope_AD nvdll_AD(i) ', my_id, nope_AD, nvdll_AD(i)
                endif
            enddo
            !write(6,*) 'my_id  nope_AD_coupl_ID ', my_id, nope_AD_coupl_ID
        else
            nope_AD_coupl_ID = nope_AD
        endif
        !write(6,*) 'my_id nope_AD_coupl_ID nvdll_AD(i) ', my_id, nope_AD_coupl_ID, nvdll_AD(nope_AD_coupl_ID)
        do i = 1, neta_AD
           ! write(6,*) 'my_id i nbdv_AD(nope_AD_coupl_ID,i)',my_id, i, nbdv_AD(nope_AD_coupl_ID,i)
        enddo
 
        ! Find flux BC data
        do i = 1, nvel_AD
            if ((lbcodei_AD(i).eq.2).or.(lbcodei_AD(i).eq.12).or.(lbcodei_AD(i).eq.22)) then
                num_flux_AD = num_flux_AD + 1
            endif
        enddo
        !write(6,*) 'num_flux_AD', num_flux_AD
        do i = 1, nvel_DG
            if ((QNIN1_DG(i).eq.0).and.((lbcodei_DG(i).eq.2).or.(lbcodei_DG(i).eq.12).or.(lbcodei_DG(i).eq.22))) then
                num_flux_DG = num_flux_DG + 1
                !write(6,*) 'my_id, i, QNIN1_DG(i) ', my_id, i, QNIN1_DG(i)
            endif
        enddo
        !write(6,*) 'num_flux_DG', num_flux_DG



        !Find maximum number of elevation and flux nodes on each side
        call mpi_allreduce(neta_AD, max_neta, 1, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
        call mpi_allreduce(num_flux_DG, max_flow_nvel, 1, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
        
        !allocate these as the max eta
        ALLOCATE (temp_ESBIN1_AD(max_neta),temp_ESBIN2_AD(max_neta))
        ALLOCATE (temp_ESBIN1_DG(max_neta),temp_ESBIN2_DG(max_neta))
        ALLOCATE (temp_QNIN1_AD(max_flow_nvel),temp_QNIN2_AD(max_flow_nvel))
        ALLOCATE (temp_QNIN1_DG(max_flow_nvel),temp_QNIN2_DG(max_flow_nvel))
        ALLOCATE (coupl_DG_id_TEMP(num_procs),coupl_AD_id_TEMP(num_procs))
        ALLOCATE (coupl_DG_id(num_procs),coupl_AD_id(num_procs))
        ALLOCATE (coupl_DG_fluxnode(num_procs), coupl_AD_Netanode(num_procs))
        ALLOCATE (coupl_DG_fluxnode_T(num_procs), coupl_AD_Netanode_T(num_procs))

        ! Dummy variables used to establish processor id of
        ! the coupling boundaries and sharing of information
        do i = 1, num_procs
            coupl_DG_id_TEMP(i) = -1
            coupl_AD_id_TEMP(i) = -1
            coupl_DG_fluxnode_T(i) = -1
            coupl_AD_Netanode_T(i) = -1
        enddo
        flow_corr_fac_TEMP = -2
        coupl_dir_TEMP = -2


        do while (adc_tprev.lt.adc_tfinal) !adc_tfinal


            
            !Increase step counter
            ntsteps = ntsteps + 1

            !Reset counters and count in time
            etmininc_count = 0
            ftmininc_count = 0
            time_a_AD = dt_AD*ntsteps + statim_AD

            !Since the data structure of the two codes is similar/identical
            !we use this trick to ensure we read from the correct files
            !and correct places in these files
            if ((nelev_bc_ad).and.(nope_AD.gt.0).and.((time_a_AD+dt_AD).gt.ETIME2_AD)) then
                etmininc_count =nint((time_a_AD+dt_AD)/etiminc_AD)*neta_AD + 1
                write(6,*) 'AD etmininc_count:  ', etmininc_count
            endif
            if ((num_flux_AD.gt.0).and.((time_a_AD+dt_AD).gt.QTIME2_AD)) then
                ftmininc_count =nint((time_a_AD+dt_AD)/ftiminc_AD)*neta_AD + 1
                write(6,*) 'AD ftmininc_count:  ', ftmininc_count
            endif
            !write(6,*) 'AD nelev_bc_ad (true/false):  ', nelev_bc_ad

            
            if ((nelev_bc_ad).and.(nope_AD.gt.0).and.(etmininc_count.ne.0)) then
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

 
            if (my_id.eq.0) then
                proc_id = num_procs
            else
                proc_id = my_id
            endif

            !Find the processors that contain the interface in each code
            ! my_id counts from ZERO!!!!
            if((ntsteps.eq.1).and.(nelev_bc_ad).and.(nope_AD.gt.0).and.(ESBIN1_AD(1).eq.0.d0)) then
                if(my_id.eq.0) then
                    coupl_AD_id_TEMP(num_procs) = my_id
                    coupl_AD_Netanode_T(num_procs) = num_neta19_AD
                else
                    coupl_AD_id_TEMP(my_id) = my_id
                    coupl_AD_Netanode_T(my_id) = num_neta19_AD
                endif
                !write(6,*) 'coupl_AD_id inside if ', coupl_AD_id_TEMP !! coupl_AD_id = 0
            endif
            if((ntsteps.eq.1).and.(num_flux_DG.gt.0).and.(QNIN1_DG(1).eq.0.d0)) then
                if(my_id.eq.0) then
                    coupl_DG_id_TEMP(num_procs) = my_id
                    coupl_DG_fluxnode_T(num_procs) = num_flux_DG
                else
                    coupl_DG_id_TEMP(my_id) = my_id
                    coupl_DG_fluxnode_T(my_id) = num_flux_DG
                endif
                !write(6,*) 'coupl_DG_id inside if ', coupl_DG_id_TEMP
            endif
            if(ntsteps.eq.1) then
                call mpi_allreduce(coupl_DG_id_TEMP, coupl_DG_id, num_procs, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
                call mpi_allreduce(coupl_AD_id_TEMP, coupl_AD_id, num_procs, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
                call  mpi_allreduce(coupl_AD_Netanode_T,coupl_AD_Netanode, num_procs, mpi_int, MPI_MAX, MPI_COMM_WORLD, ierr) ! 3 is temp
                call  mpi_allreduce(coupl_DG_fluxnode_T, coupl_DG_fluxnode, num_procs, mpi_int, MPI_MAX, MPI_COMM_WORLD, ierr) ! 1 is temp
                do i = 1, num_procs
                    

                    if (coupl_DG_fluxnode(proc_id).eq.coupl_AD_Netanode(i).and.(coupl_DG_fluxnode(proc_id).gt.0)) then
                        if (i.eq.num_procs) then
                            coupl_AD_id(proc_id) = 0
                        else
                            coupl_AD_id(proc_id) = i
                        endif
                    endif
                    if (coupl_AD_Netanode(proc_id).eq.coupl_DG_fluxnode(i).and.(coupl_AD_Netanode(proc_id).gt.0)) then
                        coupl_DG_id(proc_id) = i
                    endif
                enddo
            endif

            if(ntsteps.eq.1) then
                do i = 1, num_procs
                    if(coupl_DG_id(i).gt.0) then
                        procsum = procsum + 1
                    endif
                enddo
                write(6,*) 'Number of procs on which we couple', procsum

                if (my_id.ge.0) then
                do i = 1, num_procs
                    write(6,*) 'my_id,  coupl_DG_id', my_id, coupl_DG_id(i)
                    write(6,*) 'my_id,  coupl_DG_fluxnode', my_id,coupl_DG_fluxnode(i)
                enddo
                !if (my_id.eq.coupl_DG_id(i))
                do i = 1, num_procs
                    write(6,*) 'my_id,  coupl_AD_id', my_id, coupl_AD_id(i)
                    write(6,*) 'my_id,  coupl_AD_Netanode', my_id,coupl_AD_Netanode(i)
                    if (proc_id.eq.coupl_AD_id(i)) then
                        do j = 1,coupl_AD_Netanode(i)
                            write(6,*) 'Coupl node num', nbdv_AD(nope_AD_coupl_ID,j)
                        enddo
                    endif
                enddo
                !if (my_id.eq.coupl_AD_id(my_id))
                endif
                if (my_id.eq.0) then
                write(6,*) 'my_id,  num_flux_DG', my_id, num_flux_DG
                write(6,*) 'my_id,  neta_DG', my_id, neta_DG
                write(6,*) 'my_id,  num_flux_AD', my_id, num_flux_AD
                write(6,*) 'my_id,  neta_AD', my_id, neta_AD
                endif
            endif

            if(ntsteps.eq.1) then
                !Share interface data!
            endif
            

            ! Find normal vector on DG side to determine normal vector
            ! in the ADCIRC domain and use it to ensure
            ! correct flow direction over the coupling boundary
            if((ntsteps.eq.1).and.(proc_id.eq.coupl_DG_id(proc_id))) then
                global_edge_num = nfedn_DG(2)
                temp_NX = 0.d0
                temp_NY = 0.d0
                !write(6,*) 'global flow edge num DG ', global_edge_num
                ! no edges is num nodes -1 !
                do i = 1, num_flux_DG-1
                    !write(6,*) 'global flow edge num DG ', nfedn_DG(i)
                    !write(6,*) 'NX ', norm_X_DG(nfedn_DG(i))
                    !write(6,*) 'NY', norm_Y_DG(nfedn_DG(i))
                    temp_NY = temp_NY + norm_Y_DG(nfedn_DG(i))
                    temp_NX = temp_NX + norm_X_DG(nfedn_DG(i))
                enddo
                !write(6,*) 'temp NX ', temp_NX
                !write(6,*) 'temp NY ', temp_NY
                avg_NX = temp_NX/(num_flux_DG-1)
                avg_NY = temp_NY/(num_flux_DG-1)
                !write(6,*) 'avg NX ', avg_NX
                !write(6,*) 'avg NY ', avg_NY
                if(abs(avg_NX).gt.abs(avg_NY)) then
                    coupl_dir_TEMP = 0
                    write(6,*) 'Coupling is in the x direction, my_id ', my_id
                else if(abs(avg_NX).lt.abs(avg_NY)) then
                    coupl_dir_TEMP = 1
                    write(6,*) 'Coupling is in the y direction, my_id ', my_id
                else
                    write(6,*) 'Coupling is 45 degrees, check boundary!'
                    write(6,*)  'Calling mpi abort '
                    call mpi_abort(MPI_COMM_WORLD, 2, ierr)
                endif
                
                !Now we set the correction factor
                
                if((avg_NX.gt.0).and.(coupl_dir_TEMP.eq.0)) then
                    flow_corr_fac_TEMP = -1
                else if((avg_NX.lt.0).and.(coupl_dir_TEMP.eq.0)) then
                    flow_corr_fac_TEMP = 1
                else if((avg_NY.gt.0).and.(coupl_dir_TEMP.eq.1)) then
                    flow_corr_fac_TEMP = -1
                else if((avg_NY.lt.0).and.(coupl_dir_TEMP.eq.1)) then
                    flow_corr_fac_TEMP = 1
                endif
                write(6,*)  'Flow correction ', flow_corr_fac_TEMP
            endif
            if(ntsteps.eq.1) then
                call mpi_allreduce(flow_corr_fac_TEMP, flow_corr_fac, 1, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
                call mpi_allreduce(coupl_dir_TEMP, coupl_dir, 1, mpi_int, MPI_MAX, MPI_COMM_WORLD,ierr)
            endif
            !write(6,*) 'flow_corr_fac ', flow_corr_fac
            !write(6,*) 'coupl_dir ', coupl_dir

            !Call ADCIRC to do one step in time
            call adcirc_run(AD_ntstep)
            adc_tprev = adc_tprev + dt_AD
            !write(6,*) 'called Adcirc step', ntsteps
            !write(6,*) 'about to compute  ADCIRC outflow w nope_AD', nope_AD
            !write(6,*) 'nelev_bc_ad ,proc_id,  coupl_AD_id(proc_id) ', nelev_bc_ad, proc_id,coupl_AD_id(proc_id)
            if((nope_AD.gt.0).and.(nelev_bc_ad).and.(my_id.eq.coupl_AD_id(proc_id))) then

                !multiply veolcity by depth to get correct
                !inflow BC and the correction factor
                ! Need to ascertain correct segment in nope_AD!!!
                if (coupl_dir.eq.0) then
                    !write(6,*) 'pre ADCIRC loop nvdll_AD(nope_AD_coupl_ID)', nvdll_AD(nope_AD_coupl_ID)
                    j = 1
                    do i = 1, nvdll_AD(nope_AD_coupl_ID)
                        if ((nbfr_AD.gt.0).and.(emo_AD(1,i).eq.0)) then
                            temp_QNIN1_DG(j) = uu2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            temp_QNIN2_DG(j) = uu2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            if(eta_AD(nbdv_AD(nope_AD_coupl_ID,i)).lt.-999.d0) then
                                temp_QNIN1_DG(j) = 0.d0
                                temp_QNIN2_DG(j) = 0.d0
                            endif
                            j = j + 1
                        else if ((nbfr_AD.eq.0)) then
                            temp_QNIN1_DG(j) = uu2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            temp_QNIN2_DG(j) = uu2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            if(eta_AD(nbdv_AD(nope_AD_coupl_ID,i)).lt.-999.d0) then
                                temp_QNIN1_DG(j) = 0.d0
                                temp_QNIN2_DG(j) = 0.d0
                            endif
                            j = j + 1
                        endif
                        write(6,*) 'temp_QNIN1_DG(2)', temp_QNIN1_DG(2)
                    enddo
                else
                    j = 1
                    do i = 1, nvdll_AD(nope_AD_coupl_ID) !!Need to check correct EMO!!
                        if ((nbfr_AD.gt.0).and.(emo_AD(1,i+neta_AD-nvdll_AD(nope_AD_coupl_ID)).eq.0)) then
                            !write(6,*) 'nvdll_AD(nope_AD_coupl_ID) i j ',nvdll_AD(nope_AD_coupl_ID), i, j
                            !write(6,*) 'nbdv_AD(nope_AD_coupl_ID,i) ',nbdv_AD(nope_AD_coupl_ID,i)
                            temp_QNIN1_DG(j) = vv2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            temp_QNIN2_DG(j) = vv2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            if(eta_AD(nbdv_AD(nope_AD_coupl_ID,i)).lt.-999.d0) then
                                temp_QNIN1_DG(j) = 0.d0
                                temp_QNIN2_DG(j) = 0.d0
                            endif
                            j = j + 1
                        else if ((nbfr_AD.eq.0)) then
                            !write(6,*) 'nbdv_AD(nope_AD,i)', nbdv_AD(nope_AD,i)
                            temp_QNIN1_DG(j) = vv2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            temp_QNIN2_DG(j) = vv2_AD(nbdv_AD(nope_AD_coupl_ID,i))*flow_corr_fac*(eta_AD(nbdv_AD(nope_AD_coupl_ID,i))+dp_AD(nbdv_AD(nope_AD_coupl_ID,i)))
                            if(eta_AD(nbdv_AD(nope_AD_coupl_ID,i)).lt.-999.d0) then
                                temp_QNIN1_DG(j) = 0.d0
                                temp_QNIN2_DG(j) = 0.d0
                            endif
                            j = j + 1
                        endif
                        !write(6,*) 'temp_QNIN2_DG', temp_QNIN2_DG(i)
                    enddo
                    !write(6,*) 'pre ADCIRC loop complete i,j', i,j
                !write(6,*) 'computed ADCIRC outflow step', ntsteps
                endif
                

                !Now need to put ADCIRC data into DG code.
                !In case we are on the same processor set the
                !variables directly
                if((num_neta19_AD.eq.num_flux_DG).and.(num_flux_DG.ne.0)) then
                    j = 1
                    !write(6,*) 'pre set DG loop num_neta19_AD num_flux_DG nvel_DG',num_neta19_AD ,num_flux_DG, nvel_DG
                    do i = 1, nvel_DG
                        if ((lbcodei_DG(i).eq.2).or.(lbcodei_DG(i).eq.12).or.(lbcodei_DG(i).eq.22)) then
                            !write(6,*) 'i , j',i,j
                            QNIN2_DG(i) = temp_QNIN2_DG(j)
                            QNIN1_DG(i) = temp_QNIN2_DG(j)
                            j = j + 1
                        endif
                    enddo
                    QTIME2_DG = QTIME2_DG + dt_DG
                !write(6,*) 'shared flow to DG by setting QNIN'
                !endif
                !If we are not on the same processor we send/recieve
                !!!coupl_AD_id(proc_id) = 3
                else if((my_id.eq.coupl_AD_id(proc_id)).and.(num_neta19_AD.ne.num_flux_DG)) then
                    call mpi_send( temp_QNIN2_DG, num_neta19_AD ,mpi_real, coupl_DG_id(proc_id),10,MPI_COMM_WORLD,ierr) !coupl_DG_id = 1
                    !write(6,*) 'shared flow to DG step', ntsteps
                    !write(6,*) 'temp_QNIN1_DG(2) sent ', temp_QNIN2_DG(2)
                endif
        
            !endif
            !write(6,*) 'about to recieve flow from  ADCIRC num_flux_DG', num_flux_DG
            !write(6,*) 'num_neta19_AD , my_id,  coupl_AD_id(proc_id) ', num_neta19_AD, my_id,coupl_DG_id(proc_id)
            else if((num_flux_DG.gt.0).and.(my_id.eq.coupl_DG_id(proc_id)).and.(num_neta19_AD.eq.num_flux_DG)) then
                call mpi_recv(temp_QNIN2_DG, max_flow_nvel, mpi_real, coupl_AD_id(proc_id), 10, MPI_COMM_WORLD, status, ierr) !coupl_AD_id = 3
                !write(6,*) 'temp_QNIN1_DG(2) recieved ', temp_QNIN2_DG(2)
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

            !write(6,*) 'check outfiles step', ntsteps
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

            !Output for debugging
            !write(6,*) 'dg_dir:  ', dg_dir
            !write(6,*) 'dg_rootdir: ', trim(dg_rootdir)
            !write(6,*) 'INPUTDIR:  ', trim(dg_indir)
            !write(6,*) 'GLOBALDIR:  ', trim(dg_globdir)
            !write(6,*) 'LOCALDIR:  ', trim(dg_locdir)
            !write(6,*) 'GBLINPUTDIR:  ', trim(dg_glbindir)
            !write(6,*) 'AD LOCALDIR:  ', trim(ad_locdir)

            
            call dgswem_run
            !write(6,*) 'called DG-SWEM step', ntsteps


            !trick DG-SWEM to only step once at each call
            dg_tprev = dg_tprev + dt_DG
            iths_DG = iths_DG + 1

            if((num_flux_DG.gt.0).and.(my_id.eq.coupl_DG_id(proc_id))) then
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
                if(coupl_AD_id(proc_id).eq.coupl_DG_id(proc_id)) then
                    do i = 1, num_neta19_AD
                        ESBIN2_AD(i) = temp_ESBIN2_AD(i)
                        ESBIN1_AD(i) = temp_ESBIN1_AD(i)
                    enddo
                    ETIME2_AD = ETIME2_AD + dt_AD
                !endif
                !write(6,*) 'Trying to share elevation to ADCIRC vi hardcoding'
                else if((my_id.eq.coupl_DG_id(proc_id))) then
                    !write(6,*) 'Trying to share elevation to ADCIRC via MPI'
                    call mpi_send( temp_ESBIN2_AD, num_flux_DG ,mpi_real, coupl_AD_id(proc_id),10,MPI_COMM_WORLD,ierr)
                endif
            !write(6,*) 'shared elevation to ADCIRC, to proc coupl_AD_id(proc_id)',  coupl_AD_id(proc_id)
            !endif
            !if we are on separate processors
            else if((nope_AD.gt.0).and.(my_id.eq.coupl_AD_id(proc_id)).and.(num_neta19_AD.ne.num_flux_DG)) then
                !write(6,*) 'elevation on ADCIRC side, from proc coupl_DG_id(proc_id)', coupl_DG_id(proc_id)
                call mpi_recv(temp_ESBIN2_AD,neta_AD, mpi_real, coupl_DG_id(proc_id), 10, MPI_COMM_WORLD, status, ierr)
                !write(6,*) 'reciecved elevation from DG-SWEM '
                do i = 1, num_neta19_AD
                    ESBIN2_AD(i) = temp_ESBIN2_AD(i)
                    ESBIN1_AD(i) = temp_ESBIN2_AD(i)
                    !write(6,*) 'rec temp_ESBIN2_AD', temp_ESBIN2_AD(i))
                enddo
                ETIME2_AD = ETIME2_AD + dt_AD
                !write(6,*) 'reciecved elevation from DG-SWEM step', ntsteps
            endif
        !write(6,*) 'step  ', ntsteps

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

