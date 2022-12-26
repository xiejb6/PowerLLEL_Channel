module mod_parameters
    use mod_type
    use mod_mpi, only: myrank, ierr
    
    implicit none    
    ! make everything public unless declared private
    public

    ! Often used constants
    real(fp), parameter :: pi = acos(-1._fp)
    
    ! parallel computing parameters
    integer, save :: p_row
    integer, save :: p_col
    integer, save :: nthreads = 1

    ! mesh parameters
    logical,  save :: read_mesh
    ! mesh_type, valid only when read_mesh=.FALSE.,
    ! 0: uniform, 
    ! 1: nonuniform, clustered at the lower end, 
    ! 2: nonuniform, clustered at both ends, 
    ! 3: nonuniform, clustered at the middle
    ! stretch_ratio, valid only when mesh_type/=0, should not be zero, or the mesh becomes uniform.
    integer,  save :: mesh_type
    real(fp), save :: stretch_ratio
    integer,  save :: stretch_func
    integer,  save :: nx, ny, nz
    real(fp), save :: lx, ly, lz
    integer, parameter, dimension(6) :: nhalo = (/1, 1, 1, 1, 1, 1/)
    integer, parameter, dimension(6) :: nhalo_one = (/1, 1, 1, 1, 1, 1/)
    
    ! time step parameters
    integer,  save :: nt_end
    real(fp), save :: dt
    integer,  save :: nt_check
    real(fp), save :: cfl_limit
    real(fp), save :: div_limit

    ! restart computing parameters
    logical,       save :: is_restart
    character(60), save :: fn_prefix_input_inst = 'save'
    character(60), save :: fn_prefix_input_stat = 'save_stat'
    
    ! physical property parameters
    real(fp), save :: re, u_ref, l_ref, re_inv
    character(3), save :: initial_field
    real(fp), save :: u0 = 0.0_fp
    logical,  save :: init_with_noise = .false.
    real(fp), save :: noise_intensity = 0.0_fp
    logical,  save :: smooth_wall_visc = .false.
    real(fp), save :: u_crf = 0.0_fp

    ! forced flow parameters
    logical(1),  save, dimension(3) :: is_forced
    real(fp), save, dimension(3) :: vel_force
    
    ! statistics parameters
    integer, save :: nt_init_stat
    integer, save :: sample_interval = 0
    type(stat_info_t), save :: stat_info
    logical, save :: stat_which_var(11)

    ! data output parameters
    integer, save :: nt_out_scrn
    integer, save :: nt_out_inst
    integer, save :: nt_out_stat
    integer, save :: nt_out_save
    logical, save :: overwrite_save = .true.
    logical, save :: auto_cleanup = .true.
    integer, save :: num_retained = 0
    integer, save :: nt_out_moni
    character(4), save :: fn_prefix_inst = 'inst'
    character(4), save :: fn_prefix_save = 'save'
    character(4), save :: fn_prefix_stat = 'stat'

    ! monitor parameters
    logical, save :: out_forcing
    logical, save :: out_probe_point
    integer, save :: probe_ijk(3)
    logical, save :: out_skf_z(2)
    logical, save :: out_region = .false.
    integer, save :: region_ijk(6)
    integer, save :: nt_out_region
    character(30), save :: fn_forcing = 'monitor_forcing.out'
    character(30), save :: fn_probe = 'monitor_probe.out'
    character(30), dimension(2), save :: fn_skf_z = (/'monitor_skf_z1.out','monitor_skf_z2.out'/)
    character(30), save :: fn_prefix_region = 'monitor_region'
    
    contains
    subroutine readInputParam(fn_input)
        implicit none
        include 'mpif.h'
        character(*), intent(in) :: fn_input
        
        logical :: alive
        integer :: ios
        integer :: id_input = 10
        
        namelist /PARA/       p_row, p_col
        namelist /MESH/       read_mesh, mesh_type, stretch_ratio, stretch_func, nx, ny, nz, lx, ly, lz
        namelist /TIME/       nt_end, dt, nt_check, cfl_limit, div_limit
        namelist /RESTART/    is_restart, fn_prefix_input_inst, fn_prefix_input_stat
        namelist /PHYSICAL/   re, u_ref, l_ref, initial_field, u0, init_with_noise, noise_intensity, &
                              smooth_wall_visc, u_crf
        namelist /FORCE/      is_forced, vel_force
        namelist /STATISTICS/ nt_init_stat, sample_interval, stat_which_var
        namelist /OUTPUT/     nt_out_scrn, nt_out_inst, nt_out_stat, &
                              nt_out_save, overwrite_save, auto_cleanup, num_retained, &
                              nt_out_moni

        inquire(file=fn_input, exist=alive)
        if (.not. alive) then
            if (myrank == 0) then
                write(*,'(A)') "PowerLLEL.ERROR.readInputParam: File param.in doesn't exist!"
            endif
            call MPI_FINALIZE(ierr)
            stop
        endif
        open(unit=id_input, file=fn_input, status='old', iostat=ios)
            if (ios /= 0) then
                if (myrank == 0) then
                    write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Fail to open file param.in!"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif

            read(id_input, nml=PARA,       end=20, err=21, iostat=ios)
            read(id_input, nml=MESH,       end=20, err=22, iostat=ios)
            read(id_input, nml=TIME,       end=20, err=23, iostat=ios)
            read(id_input, nml=RESTART,    end=20, err=24, iostat=ios)
            read(id_input, nml=PHYSICAL,   end=20, err=25, iostat=ios)
            read(id_input, nml=FORCE,      end=20, err=30, iostat=ios)
            read(id_input, nml=STATISTICS, end=20, err=31, iostat=ios)
            read(id_input, nml=OUTPUT,     end=20, err=32, iostat=ios)

    20      rewind(id_input)
    21      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with PARA line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
    22      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with MESH line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
    23      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with TIME line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
    24      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with RESTART line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
    25      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with PHYSICAL line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
    30      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with FORCE line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
    31      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with STATISTICS line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
    32      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readInputParam: Problem with OUTPUT line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
        close(id_input)

        ! read parameters of the monitor
        call readMonitorParam(fn_input)

        ! do a sanity check on input parameters
        call checkInputParam()

        re_inv = u_ref*l_ref/re
        
        return
    end subroutine readInputParam

    subroutine readMonitorParam(fn_input)
        implicit none
        include 'mpif.h'
        character(*), intent(in) :: fn_input

        integer :: ig, jg, kg
        logical :: alive
        integer :: ios
        integer :: id_input = 10

        namelist /MONITOR/ out_forcing, out_probe_point, probe_ijk, out_skf_z, &
                           out_region, nt_out_region, region_ijk

        inquire(file=fn_input, exist=alive)
        if (.not. alive) then
            if (myrank == 0) then
                write(*,'(A)') "PowerLLEL.ERROR.readMonitorParam: File "//trim(adjustl(fn_input))//" doesn't exist!"
            endif
            call MPI_FINALIZE(ierr)
            stop
        endif
        open(unit=id_input, file=fn_input, status='old', iostat=ios)
            if (ios /= 0) then
                if (myrank == 0) then
                    write(*,'(A)') "PowerLLEL.ERROR.readMonitorParam: Fail to open file "//trim(adjustl(fn_input))//"!"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif

            read(id_input, nml=MONITOR, end=20, err=21, iostat=ios)

    20      rewind(id_input)
    21      if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.readMonitorParam: Problem with MONITOR line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
        close(id_input)

        return
    end subroutine readMonitorParam

    subroutine checkInputParam()
        implicit none
        include 'mpif.h'

        logical :: check_passed
        logical :: subcheck_passed

        check_passed = .true.

        call checkMeshAndPara(subcheck_passed)
        check_passed = check_passed .and. subcheck_passed

        if (.not. check_passed) then
            call MPI_FINALIZE(ierr)
            stop
        endif

        return
    end subroutine checkInputParam

    subroutine checkMeshAndPara(subcheck_passed)

        implicit none

        logical, intent(inout) :: subcheck_passed
        logical :: passed

        subcheck_passed = .true.

        passed = ( mod(nx,2)==0 .and. mod(ny,2)==0 .and. mod(nz,2)==0 )
        if ( (.not.passed) .and. myrank==0 ) &
            write(*,'(A)') "PowerLLEL.ERROR.checkInputParam: <nx>, <ny>, <nz> should be even integers!"
        subcheck_passed = subcheck_passed .and. passed

        passed = ( mod(nx,p_row)==0 )
        if ( (.not.passed) .and. myrank==0 ) &
            write(*,'(A)') "PowerLLEL.ERROR.checkInputParam: <nx> should be exactly divisible by <p_row>!"
        subcheck_passed = subcheck_passed .and. passed

#ifdef _PDD
        passed = ( mod(ny,p_row)==0)
        if ( (.not.passed) .and. myrank==0 ) &
            write(*,'(A)') "PowerLLEL.ERROR.checkInputParam: <ny> should be exactly divisible by <p_row>!"
#else
        passed = ( mod(ny,p_row)==0 .and. mod(ny,p_col)==0 )
        if ( (.not.passed) .and. myrank==0 ) &
            write(*,'(A)') "PowerLLEL.ERROR.checkInputParam: <ny> should be exactly divisible by both <p_row> and <p_col>!"
#endif
        subcheck_passed = subcheck_passed .and. passed

        passed = ( mod(nz,p_col)==0 )
        if ( (.not.passed) .and. myrank==0 ) &
            write(*,'(A)') "PowerLLEL.ERROR.checkInputParam: <nz> should be exactly divisible by <p_col>!"
        subcheck_passed = subcheck_passed .and. passed
        
        passed = ( mod(nz,p_col*nthreads)==0 )
        if ( (.not.passed) .and. myrank==0 ) &
            write(*,'(A)') "PowerLLEL.ERROR.checkInputParam: <nz> should be exactly divisible by <p_col*nthreads>!"
        subcheck_passed = subcheck_passed .and. passed
        
        return
    end subroutine checkMeshAndPara
    
end module mod_parameters