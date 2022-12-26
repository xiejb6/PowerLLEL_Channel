program PowerLLEL_Convert_Precision
    use mod_dataio_postproc, only: inputField, outputField
    use mod_hdf5, only: initIO, freeIO
    use mod_mpi, only: st, sz, myrank, ierr, initMPI, freeMPI

    implicit none

    include 'mpif.h'

#ifndef R4_TO_R8
    integer, parameter :: fp_src = selected_real_kind(15)
    integer, parameter :: fp_dst = selected_real_kind(6)
#else
    integer, parameter :: fp_src = selected_real_kind(6)
    integer, parameter :: fp_dst = selected_real_kind(15)
#endif
    integer, parameter, dimension(6) :: nhalo = (/0, 0, 0, 0, 0, 0/)
    integer :: nx, ny, nz, p_row, p_col
    integer :: nt = -1
    character(80) :: str, prefix, var

    ! main arrays
    real(fp_src), allocatable, dimension(:,:,:), save :: var_src
    real(fp_dst), allocatable, dimension(:,:,:), save :: var_dst

    double precision :: wtime, tmp
    integer :: nranks
    integer :: istatus
    integer :: i, j, k
    character :: byte
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    if (myrank == 0) then
        write(*,'(A)') "================================================================================"
        write(*,'(A)') "                           PowerLLEL_Convert_Precision                          "
        write(*,'(A)') "================================================================================"
        write(*,'(A)') "PowerLLEL_Convert_Precision.NOTE: Start postprocessing ..."
    endif

    ! read input parameters from command line
    CALL get_command_argument(1, str)
    if (len_trim(str) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Convert_Precision.NOTE: Please pass grid number <nx> as an argument to PowerLLEL_Convert_Precision!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    read(str, '(I5)') nx
    CALL get_command_argument(2, str)
    if (len_trim(str) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Convert_Precision.NOTE: Please pass grid number <ny> as an argument to PowerLLEL_Convert_Precision!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    read(str, '(I5)') ny
    CALL get_command_argument(3, str)
    if (len_trim(str) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Convert_Precision.NOTE: Please pass grid number <nz> as an argument to PowerLLEL_Convert_Precision!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    read(str, '(I5)') nz
    CALL get_command_argument(4, str)
    if (len_trim(str) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Convert_Precision.NOTE: Please pass grid number <p_row> as an argument to PowerLLEL_Convert_Precision!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    read(str, '(I5)') p_row
    CALL get_command_argument(5, str)
    if (len_trim(str) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Convert_Precision.NOTE: Please pass grid number <p_col> as an argument to PowerLLEL_Convert_Precision!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    read(str, '(I5)') p_col

    CALL get_command_argument(6, prefix)
    if (len_trim(prefix) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Convert_Precision.NOTE: Please pass an inst. file prefix as an argument to PowerLLEL_Convert_Precision!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    CALL get_command_argument(7, var)
    if (len_trim(var) == 0) then
        if (myrank == 0) write(*,'(A)') &
        "PowerLLEL_Convert_Precision.NOTE: Please pass an inst. file variable (e.g., u) as an argument to PowerLLEL_Convert_Precision!"
        call MPI_FINALIZE(ierr)
        stop
    endif
    if (myrank == 0) write(*,'(A)') "PowerLLEL_Convert_Precision.NOTE: Inst. file is '"//trim(prefix)//trim(var)//".h5'"


    ! initialize MPI
    call initMPI(nx, ny, nz, (/'PP','PP','NN'/), p_row, p_col, nhalo)

    ! initialize parallel IO
    call initIO(MPI_COMM_WORLD)

    ! allocate variables
    allocate(var_src(1-nhalo(1):sz(1)+nhalo(2), 1-nhalo(3):sz(2)+nhalo(4), 1-nhalo(5):sz(3)+nhalo(6)), stat=istatus)
    if (istatus /= 0) call abort(101, "PowerLLEL_Convert_Precision: Out of memory!")
    allocate(var_dst(1-nhalo(1):sz(1)+nhalo(2), 1-nhalo(3):sz(2)+nhalo(4), 1-nhalo(5):sz(3)+nhalo(6)), stat=istatus)
    if (istatus /= 0) call abort(101, "PowerLLEL_Convert_Precision: Out of memory!")

    ! read fields
    call inputField(trim(prefix)//trim(var)//'.h5', nt, st, sz, trim(var), var_src)

    ! convert precision
    if (myrank == 0) write(*,'(A)') "PowerLLEL_Convert_Precision.NOTE: Start converting ..."
    wtime = MPI_WTIME()

    do k = 1, sz(3)
        do j = 1, sz(2)
            do i = 1, sz(1)
                var_dst(i, j, k) = real(var_src(i, j, k), fp_dst)
            enddo
        enddo
    enddo

    wtime = MPI_WTIME() - wtime
    call MPI_REDUCE(wtime, tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks, ierr)
    wtime = tmp/nranks
    if (myrank == 0) write(*,'(A,1PE10.3,A)') &
    "PowerLLEL_Convert_Precision.NOTE: Finish converting in ", wtime, "s"

    ! write fields
    write(byte, '(I1)') fp_dst
    call outputField('real'//byte//'_'//trim(prefix), &
                     nt, (/nx,ny,nz/), st, sz, nhalo, trim(var), var_dst)

    deallocate(var_src)
    deallocate(var_dst)

    call freeIO()
    call freeMPI()

    if (myrank == 0) then
        write(*,'(A)') "********************************************************************************"
        write(*,'(A)') "PowerLLEL_Convert_Precision.NOTE: Postprocessing ends successfully!"
    endif

    call MPI_FINALIZE(ierr)
    
end program PowerLLEL_Convert_Precision