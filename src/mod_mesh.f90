module mod_mesh
    use mod_type,       only: fp
    use mod_parameters, only: pi, read_mesh, mesh_type, stretch_ratio, stretch_func, lx, ly, lz, &
                              nx_global=>nx, ny_global=>ny, nz_global=>nz, smooth_wall_visc
    use mod_mpi,        only: myrank, ierr, sz, st, neighbor
    use mod_hdf5,       only: HID_T, createFile, write1d, closeFile

    implicit none
    ! make everything public unless declared private
    public

    ! arrays for mesh
    ! 'c' refers to '(mesh cell) center', 'f' refers to '(mesh cell) face'
    real(fp), save  :: dx, dy
    real(fp), save  :: dx_inv, dy_inv
    real(fp), allocatable, dimension(:), save :: xc_global   ! redundant, but used in the output process
    real(fp), allocatable, dimension(:), save :: yc_global   ! redundant, but used in the output process
    real(fp), allocatable, dimension(:), save :: zc_global   ! redundant, but used in the output process
    real(fp), allocatable, dimension(:), save :: dzf_global
    real(fp), allocatable, dimension(:), save :: dzc, dzc_inv
    real(fp), allocatable, dimension(:), save :: dzf, dzf_inv
    real(fp), allocatable, dimension(:), save :: dzflzi
    real(fp), allocatable, dimension(:), save :: visc_dzf_inv

    contains
    subroutine initMesh()
        implicit none
        include 'mpif.h'

        integer :: i, j, k
        
        allocate(xc_global(0:nx_global+1))
        allocate(yc_global(0:ny_global+1))
        allocate(zc_global(0:nz_global+1))
        allocate(dzf_global(0:nz_global+1))
        allocate(dzc(0:sz(3)+1), dzc_inv(0:sz(3)+1))
        allocate(dzf(0:sz(3)+1), dzf_inv(0:sz(3)+1))
        allocate(dzflzi(0:sz(3)+1))
        allocate(visc_dzf_inv(0:sz(3)+1))

        ! mesh in x & y direction
        dx = lx/nx_global
        dy = ly/ny_global
        dx_inv = 1.0_fp/dx
        dy_inv = 1.0_fp/dy
        xc_global(0) = -0.5_fp*dx
        do i = 1, nx_global+1; xc_global(i) = xc_global(0) + i*dx; enddo
        yc_global(0) = -0.5_fp*dy
        do j = 1, ny_global+1; yc_global(j) = yc_global(0) + j*dy; enddo

        ! mesh in z direction
        if (read_mesh) then
            ! read the mesh spacing <dzf_global>
            call inputMesh(nz_global, dzf_global)
        else
            ! generate <dzf_global> by mesh functions
            call initMeshByFunc(mesh_type, stretch_ratio, nz_global, lz, dzf_global)
        endif
        
        zc_global(0) = -0.5_fp*dzf_global(0)
        do k = 1, nz_global+1
            zc_global(k) = zc_global(k-1) + 0.5_fp*(dzf_global(k-1) + dzf_global(k))
        enddo
        do k = 0, sz(3)+1
            dzf(k) = dzf_global(st(3)-1+k)
            dzf_inv(k) = 1.0_fp/dzf(k)
        enddo
        do k = 0, sz(3)
            dzc(k) = 0.5_fp*(dzf(k)+dzf(k+1))
            dzc_inv(k) = 1.0_fp/dzc(k)
        enddo
        dzflzi(:) = dzf(:)/lz

        visc_dzf_inv(:) = dzf_inv(:)
        if (smooth_wall_visc) then
            if (neighbor(5) == MPI_PROC_NULL) visc_dzf_inv(1) = 0.5_fp * (dzc(0) + dzc(1))
            if (neighbor(6) == MPI_PROC_NULL) visc_dzf_inv(sz(3)) = 0.5_fp * (dzc(sz(3)-1) + dzc(sz(3)))
        endif

        ! output mesh files for post-processing and check
        call outputMesh()
        
        return
    end subroutine initMesh

    subroutine inputMesh(nz_global, dzf_global)
        implicit none
        include 'mpif.h'
        integer, intent(in) :: nz_global
        real(fp), dimension(0:), intent(out) :: dzf_global

        logical :: alive
        integer :: k, ios
        integer :: id_mesh = 11
        character(7) :: fn_mesh = 'mesh.in'

        inquire(file=fn_mesh, exist=alive)
        if (.not. alive) then
            if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.inputMesh: File "//fn_mesh//" doesn't exist!"
            call MPI_FINALIZE(ierr)
            stop
        endif
        open(unit=id_mesh, file=fn_mesh, status='old', iostat=ios)
            if (ios /= 0) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.inputMesh: Fail to open file "//fn_mesh//"!"
                call MPI_FINALIZE(ierr)
                stop
            endif
            do k = 1, nz_global
                read(id_mesh,*) dzf_global(k)
            enddo
        close(id_mesh)

        dzf_global(0) = dzf_global(1)
        dzf_global(nz_global+1) = dzf_global(nz_global)

        return
    end subroutine inputMesh

    subroutine initMeshByFunc(mesh_type, stretch_ratio, nz_global, lz, dzf_global)
        implicit none
        integer, intent(in) :: mesh_type, nz_global
        real(fp), intent(in) :: stretch_ratio, lz
        real(fp), dimension(0:), intent(out) :: dzf_global

        real(fp), dimension(0:nz_global) :: zf
        real(fp) :: z0
        integer  :: k

        if (stretch_ratio == 0.0_fp) then
            dzf_global(:) = lz/nz_global
            return
        else
            select case(mesh_type)
            case(0)
                ! uniform
                dzf_global(:) = lz/nz_global
            case(1)
                ! nonuniform, clustered at both ends
                zf(0) = 0.0_fp
                select case(stretch_func)
                case(0)
                    do k = 1, nz_global
                        z0 = k/(1.0_fp*nz_global)
                        zf(k) = 0.5_fp*( 1.0_fp + tanh((z0-0.5_fp)*stretch_ratio)/tanh(stretch_ratio/2.0_fp) )
                    enddo
                case(1)
                    do k = 1, nz_global
                        z0 = k/(1.0_fp*nz_global)
                        zf(k) = 0.5_fp*( 1.0_fp + sin((z0-0.5_fp)*stretch_ratio*pi)/sin(stretch_ratio*pi/2.0_fp) )
                    enddo
                end select
                do k = 1, nz_global
                    dzf_global(k) = lz*(zf(k) - zf(k-1))
                enddo
                dzf_global(0) = dzf_global(1)
                dzf_global(nz_global+1) = dzf_global(nz_global)
            end select
        endif

        return
    end subroutine initMeshByFunc

    subroutine outputMesh()
        implicit none
        integer(HID_T) :: fh
        character(5 ) :: string_dump
        character(30) :: fn_mesh
        integer :: id_mesh = 12
        integer :: k

        write(string_dump, '(I5)') nx_global
        fn_mesh = 'mesh_'//trim(adjustl(string_dump))//'-'
        write(string_dump, '(I5)') ny_global
        fn_mesh = trim(adjustl(fn_mesh))//trim(adjustl(string_dump))//'-'
        write(string_dump, '(I5)') nz_global
        fn_mesh = trim(adjustl(fn_mesh))//trim(adjustl(string_dump))//'.h5'
        call createFile(trim(adjustl(fn_mesh)), fh)
        call write1d(fh, nx_global, 'xc', xc_global(1:nx_global))
        call write1d(fh, ny_global, 'yc', yc_global(1:ny_global))
        call write1d(fh, nz_global, 'zc', zc_global(1:nz_global))
        call closeFile(fh)
        if (myrank == 0) &
            write(*,'(A)') "PowerLLEL.NOTE.outputMesh: Finish writing file "//trim(adjustl(fn_mesh))//"!"

        if (myrank == 0) then
            write(string_dump, '(I5.5)') nz_global
            fn_mesh = 'mesh_'//string_dump//'.out'
            open(unit=id_mesh, file=fn_mesh, status='replace')
            write(id_mesh, '(4X,A,1X,20X,A,1X,21X,A)') 'k', 'dzf', 'zc'
            do k = 1, nz_global
                write(id_mesh, '(I5.5,2(1X,E23.15))') k, dzf_global(k), zc_global(k)
            enddo
            close(id_mesh)
            write(*,'(A)') "PowerLLEL.NOTE.outputMesh: Finish writing file "//trim(adjustl(fn_mesh))//"!"
        endif
        
        return
    end subroutine outputMesh

    subroutine freeMesh()
        implicit none

        deallocate(xc_global)
        deallocate(yc_global)
        deallocate(zc_global)
        deallocate(dzf_global)
        deallocate(dzc, dzc_inv)
        deallocate(dzf, dzf_inv)
        deallocate(dzflzi)
        deallocate(visc_dzf_inv)

        return
    end subroutine freeMesh

end module mod_mesh