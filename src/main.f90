program PowerLLEL
    use, intrinsic :: iso_c_binding, only: c_null_char
    use mod_type,          only: fp
    use mod_parameters
    use mod_mpi,           only: MPI_REAL_FP, myrank, ierr, sz, initMPI, freeMPI
    use mod_mpi,           only: comm_cart, halotype_vel, halotype_one, neighbor, neighbor_xyz
    use mod_variables
    use mod_mesh
    use mod_poissonSolver, only: initPoissonSolver, executePoissonSolver, freePoissonSolver
#ifdef USE_C
    use mod_poissonSolver, only: init_poisson_solver, execute_poisson_solver, free_poisson_solver
#endif
    use mod_initFlow,      only: initFlow
    use mod_updateBound,   only: updateBoundVel, updateBoundP
    use mod_calcVel,       only: timeIntVelRK1, timeIntVelRK2, correctVel, forceVel, transform2CRF
    use mod_calcRHS,       only: calcRHS
    use mod_monitor,       only: initMonitor, freeMonitor
    use mod_hdf5,          only: initIO, freeIO
    use mod_dataIO,        only: inputData, outputData, inputStatData
    use mod_statistics,    only: allocStat, freeStat, initStat, calcStat
    use mod_utils,         only: value_index_pair_t, initCheckCFLAndDiv, freeCheckCFLAndDiv, &
                                 calcMaxCFL, calcMaxDiv, checkCFL, checkDiv, checkNaN
    !$ use omp_lib
#ifdef GPTL
    use gptl
#endif

    implicit none

    include 'mpif.h'

    integer :: nt, nt_in, nt_start
    real(fp):: wtime, wtime_avg
    type(value_index_pair_t):: cfl_max, div_max
    logical :: check_passed, subcheck_passed
    character(13) :: remaining_time_c

    integer :: ret
#ifdef GPTL

    ! GPTL initialization
    call gptlprocess_namelist('gptlnl', 1, ret)
    if (ret /= 0) call abort(999, "main: GPTL namelist read failure!")
    ret = gptlinitialize()
#endif


    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, ret, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    !$ nthreads = OMP_GET_MAX_THREADS()

    if (myrank == 0) then
        write(*,'(A)') "================================================================================"
        write(*,'(A)') "                      PowerLLEL_channel (DNS, 2ndFD + RK2)                      "
        write(*,'(A)') "================================================================================"
        write(*,'(A)') "PowerLLEL.NOTE: Initialization starts ..."
    endif

    ! read input parameters from file
    call readInputParam('param.in')

    ! initialize MPI
    call initMPI(nx, ny, nz, (/'PP','PP','NN'/), p_row, p_col, nhalo)

    ! initialize parallel IO
    call initIO(MPI_COMM_WORLD)
    
    ! initialize the monitoring point
    call initMonitor()
    
    ! allocate variables
    call allocVariables(nhalo, sz)

    ! initialize mesh
    call initMesh()

    ! initialize Poisson solver
#ifdef USE_C
    call init_poisson_solver(nx, ny, nz, dx, dy, dzf_global, "PP"//c_null_char, "PP"//c_null_char, "NN"//c_null_char, &
                             neighbor_xyz)
#else
    call initPoissonSolver(nx, ny, nz, dx, dy, dzf_global, (/'PP','PP','NN'/), real((/0.0,0.0, 0.0,0.0, 0.0,0.0/),fp))
#endif

    ! initialize the flow field
    if (is_restart) then
        if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE: Initializing flow from checkpoint fields ..."
        call inputData(nt_in, u, v, w)
        nt_start = nt_in + 1
    else
        if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE: Initializing flow according to input parameters ..."
        call initFlow(u, v, w)
        nt_start = 1
    endif
    
    ! important! subtract the convecting reference frame velocity u_crf from u
    call transform2CRF(u_crf, nhalo, sz, u, vel_force(1))

    ! update velocity & pressure boundary conditions
    call updateBoundVel(u, v, w, u_crf)
    call updateBoundP(p)

    ! load the statistics data if necessary
    if (is_restart .and. nt_start > nt_init_stat+1) then
        call allocStat(sz)
        call inputStatData(nt_start-1)
    endif

    ! initialize MPI variables related to the calculation of CFL and Divergence
    call initCheckCFLAndDiv(cfl_max, div_max)

    if (myrank == 0) then
        write(*,'(A)') "PowerLLEL.NOTE: Initialization ends successfully!"
        write(*,'(A,I9,A)') "PowerLLEL.NOTE: Simulation starts at nt = ", nt_start, "!"
        write(*,'(A)') "********************************************************************************"
        write(*,999) 'nt', 'speed(wSteps/Day)', 'remaining time', 'cfl_max', 'div_max'
    999 format(A9,2X,A17,2X,A14,2X,A10,2X,A10)
    endif

#ifdef GPTL
    ret = gptlstart('Main loop')
#endif

    !Start timing
    wtime = MPI_WTIME()

    !===========================!
    !  Main time marching loop  !
    !===========================!
    do nt = nt_start, nt_end

#ifdef GPTL
        ret = gptlstart('uvw1')
#endif

        call timeIntVelRK1(u, v, w, u1, v1, w1, u_crf)

#ifdef GPTL
        ret = gptlstop('uvw1')
        ret = gptlstart('Update boundary vel')
#endif

#ifndef NB_HALO
        call updateBoundVel(u1, v1, w1, u_crf)
#endif

#ifdef GPTL
        ret = gptlstop('Update boundary vel')
        ret = gptlstart('uvw2')
#endif

        call timeIntVelRK2(u, v, w, u1, v1, w1, u_crf)

#ifdef GPTL
        ret = gptlstop('uvw2')
        ret = gptlstart('Update boundary vel')
#endif

#ifndef NB_HALO
        call updateBoundVel(u, v, w, u_crf)
#endif

#ifdef GPTL
        ret = gptlstop('Update boundary vel')
        ret = gptlstart('Calculate RHS')
#endif

        call calcRHS(u, v, w, p)

#ifdef GPTL
        ret = gptlstop('Calculate RHS')
        ret = gptlstart('Poisson solver')
#endif

#ifdef USE_C
        call execute_poisson_solver(p)
#else
        call executePoissonSolver(p)
#endif

#ifdef GPTL
        ret = gptlstop('Poisson solver')
        ret = gptlstart('Update boundary pres')
#endif

        call updateBoundP(p)

#ifdef GPTL
        ret = gptlstop('Update boundary pres')
        ret = gptlstart('Correct vel')
#endif

        call correctVel(p, u, v, w)

#ifdef GPTL
        ret = gptlstop('Correct vel')
        ret = gptlstart('Force vel')
#endif

        call forceVel(u, v, w)

#ifdef GPTL
        ret = gptlstop('Force vel')
        ret = gptlstart('Update boundary vel')
#endif

        call updateBoundVel(u, v, w, u_crf)

#ifdef GPTL
        ret = gptlstop('Update boundary vel')
#endif

        if (mod(nt, nt_check) == 0) then
            call checkNaN(u, 'u', subcheck_passed)
            if (myrank == 0) check_passed = subcheck_passed
            call calcMaxCFL(nhalo, sz, dt, dx_inv, dy_inv, dzf_inv, u, v, w, cfl_max)
            call checkCFL(cfl_limit, cfl_max, subcheck_passed)
            if (myrank == 0) check_passed = check_passed .and. subcheck_passed
            call calcMaxDiv(nhalo, sz, dx_inv, dy_inv, dzf_inv, u, v, w, div_max)
            call checkDiv(div_limit, div_max, subcheck_passed)
            if (myrank == 0) check_passed = check_passed .and. subcheck_passed
            call MPI_BCAST(check_passed, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            if (.not. check_passed) then
                if (myrank == 0) then
                    write(*,997) 'PowerLLEL.ERROR: nt = ',nt,', cfl_max = ',cfl_max%value, &
                                 ' at (',cfl_max%ig,',',cfl_max%jg,',',cfl_max%kg,') from rank', cfl_max%rank
                    write(*,997) 'PowerLLEL.ERROR: nt = ',nt,', div_max = ',div_max%value, &
                                 ' at (',div_max%ig,',',div_max%jg,',',div_max%kg,') from rank', div_max%rank
                997 format(A,I9,A,E10.3,4(A,I5))
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        endif

        if (nt>nt_init_stat) then
            if (nt == nt_init_stat+1) then
                call allocStat(sz)
                call initStat()
                if (myrank == 0) write(*,'(A,I9,A)') "PowerLLEL.NOTE: Statistical process starts at nt = ", nt, "!"
            endif
            call calcStat(nt, u, v, w, p, u_crf)
        endif

        if (mod(nt, nt_out_scrn) == 0) then
            wtime = MPI_WTIME() - wtime
            wtime_avg = wtime
            call MPI_REDUCE(wtime, wtime_avg, 1, MPI_REAL_FP, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            wtime_avg = wtime_avg/p_row/p_col
            if (myrank == 0) then
                call convertTime(int((nt_end-nt)*wtime_avg/nt_out_scrn), remaining_time_c)
                write(*,998) nt, nt_out_scrn*3600.0*24.0/wtime_avg/10000, remaining_time_c, cfl_max%value, div_max%value
            998 format(I9,2X,F17.3,2X,A14,2(2X,E10.3))
            endif
        endif
     
        ! output routines below
        call outputData(nt, u, v, w, p, u_crf)

        if (mod(nt, nt_out_scrn) == 0) then
            wtime = MPI_WTIME()
        endif

    enddo

#ifdef GPTL
    ret = gptlstop('Main loop')
#endif

    call freeIO()
#ifdef USE_C
    call free_poisson_solver()
#else
    call freePoissonSolver()
#endif
    call freeMesh()
    call freeMonitor()
    call freeVariables()
    call freeStat()
    call freeCheckCFLAndDiv()
    call freeMPI()

#ifdef GPTL
    ! Print timer stats to file named 'timing.summary' and 'timing.$(myrank)'
    ret = gptlpr_summary(MPI_COMM_WORLD)
    ! ret = gptlpr(myrank)
    ret = gptlfinalize()
#endif

    if (myrank == 0) then
        write(*,'(A)') "********************************************************************************"
        write(*,'(A)') "PowerLLEL.NOTE: Simulation ends successfully!"
    endif

    call MPI_FINALIZE(ierr)

    contains
    subroutine convertTime(t, tc)
        implicit none
        integer, intent(in) :: t
        character(13), intent(out) :: tc

        integer :: t_d, t_h, t_m, t_s, t_res

        t_d = int(t/86400)
        t_res = t-t_d*86400
        t_h = int(t_res/3600)
        t_res = t_res-t_h*3600
        t_m = int(t_res/60)
        t_s = t_res-t_m*60

        write(tc,'(I3,A,3(I2,A))') t_d,'d',t_h,'h',t_m,'m',t_s,'s'

    end subroutine convertTime

end program PowerLLEL