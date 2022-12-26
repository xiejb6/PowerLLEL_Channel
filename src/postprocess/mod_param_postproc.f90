module mod_param_postproc
    use mod_type, only: fp
    use mod_mpi,  only: myrank, ierr
    
    implicit none
    ! make everything public unless declared private
    public

    ! Often used constants
    real(fp), parameter :: pi = acos(-1._fp)

    ! parallel computing parameters
    ! character(4), save :: decomp_mode = 'xpen'     ! 'xpen', 'ypen' and 'zpen' are available
    integer, save :: p_row
    integer, save :: p_col

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
    integer, parameter, dimension(6) :: nhalo_zero = (/0, 0, 0, 0, 0, 0/)

    ! postprocessing parameters
    logical, save :: is_vel_staggered = .true.
    logical, save :: out_vorticity = .false.
    logical, save :: out_q = .false.
    logical, save :: out_lambda2 = .false.

    contains
    subroutine readInputParam(fn_input)
        implicit none
        include 'mpif.h'
        character(*), intent(in) :: fn_input

        logical :: alive
        integer :: ios
        integer :: id_input = 10

        namelist /PARA/ &
            p_row, p_col
        namelist /MESH/ &
            read_mesh, mesh_type, stretch_ratio, stretch_func, nx, ny, nz, lx, ly, lz
        namelist /POST_VORTEX/ & 
            is_vel_staggered, out_vorticity, out_q, out_lambda2
        
        inquire(file=fn_input, exist=alive)
        if (.not. alive) then
            if (myrank == 0) then
                write(*,'(A)') "PowerLLEL_Postprocess.ERROR.readInputParam: File "//trim(fn_input)//" doesn't exist!"
            endif
            call MPI_FINALIZE(ierr)
            stop
        endif
        open(unit=id_input, file=fn_input, status='old', iostat=ios)
            if (ios /= 0) then
                if (myrank == 0) then
                    write(*,'(A)') "PowerLLEL_Postprocess.ERROR.readInputParam: Fail to open file "//trim(fn_input)//"!"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif

            read(id_input, nml=PARA,      end=20, err=21, iostat=ios)
            read(id_input, nml=MESH,      end=20, err=22, iostat=ios)
            read(id_input, nml=POST_VORTEX, end=20, err=23, iostat=ios)
        20  rewind(id_input)
        21  if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL_Postprocess.ERROR.readInputParam: Problem with PARA line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
        22  if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL_Postprocess.ERROR.readInputParam: Problem with MESH line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
        23  if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL_Postprocess.ERROR.readInputParam: Problem with POST_VORTEX line in the input file!"
                call MPI_FINALIZE(ierr)
                stop
            endif
        close(id_input)

    end subroutine readInputParam

end module mod_param_postproc