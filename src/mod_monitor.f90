module mod_monitor
    use mod_type,       only: fp
    use mod_parameters, only: out_forcing, out_probe_point, probe_ijk, out_skf_z, &
                              out_region, nt_out_region, region_ijk, &
                              fn_forcing, fn_probe, fn_skf_z, fn_prefix_region, &
                              nx, ny, nz, re_inv, u_ref, is_forced, is_restart, nhalo
    use mod_mesh,       only: dzf, dzflzi
    use mod_mpi,        only: MPI_REAL_FP, myrank, ierr, st, en, sz
    use mod_utils,      only: Mean_h, abort
    use mod_hdf5,       only: HID_T, createFile, writeAttribute, write3d, closeFile
    
    implicit none

    include 'mpif.h'
    
    ! make everything private unless declared public
    private

    integer, save :: id_forcing = 11

    type probe_point_t
        logical  :: have_this_point    ! decide whether the local mpi rank has the probe point
        integer  :: i, j, k            ! probe point local index
        integer  :: ig, jg, kg         ! probe point global index
        real(fp) :: u, v, w, p         ! physical quantities at the point
    end type
    type(probe_point_t), save :: probe_point
    integer, save :: id_probe = 12

    integer, dimension(2), save :: id_skf_z = (/13, 14/)
    integer, dimension(2), save :: comm_skf_z, myrank_skf_z

    type region_t
        logical :: have_this_region
        integer, dimension(3) :: gsize, offset
        integer, dimension(3) :: stg, eng
        integer, dimension(3) :: st, sz
    end type
    type(region_t), save :: ext_region
    real(fp), allocatable, dimension(:,:,:), save :: ext_buffer

    public :: probe_point, initMonitor, outputMonitor, freeMonitor

contains

    subroutine initMonitor()
        implicit none
        integer :: ig, jg, kg
        integer :: istatus

        if (out_forcing) then
            if (is_restart) then
                open(unit=id_forcing, file=fn_forcing, position='append')
            else
                open(unit=id_forcing, file=fn_forcing)
            endif
            if (myrank == 0) write(id_forcing,'(A9,3(A15))') 'nt', 'u_mean', 'v_mean', 'w_mean'
        endif

        if (out_probe_point) then
            ig = probe_ijk(1)
            jg = probe_ijk(2)
            kg = probe_ijk(3)

            if ( (ig >= st(1) .and. ig <= en(1)) .and. &
                (jg >= st(2) .and. jg <= en(2)) .and. &
                (kg >= st(3) .and. kg <= en(3)) ) then
                probe_point%have_this_point = .true.
                probe_point%ig = ig
                probe_point%jg = jg
                probe_point%kg = kg
                probe_point%i  = ig-st(1)+1
                probe_point%j  = jg-st(2)+1
                probe_point%k  = kg-st(3)+1
            endif
            
            if (is_restart) then
                open(unit=id_probe, file=fn_probe, position='append')
            else
                open(unit=id_probe, file=fn_probe)
            endif
            if (probe_point%have_this_point) then
                write(id_probe,'(A,3(I4,A))') 'Probe point (', probe_point%ig, ',', probe_point%jg, ',', probe_point%kg, '): '
                write(id_probe,'(A9,4(A15))') 'nt', 'u', 'v', 'w', 'p'
            endif
        endif

        ! initialize new MPI communicators for the calculation of skin friction coefficients, 
        ! at the bottom/top wall of the computational domain
        if (out_skf_z(1)) then
            call initCommForSkf( 1, st(3), en(3), comm_skf_z(1), myrank_skf_z(1))
            if (is_restart) then
                open(unit=id_skf_z(1), file=fn_skf_z(1), position='append')
            else
                open(unit=id_skf_z(1), file=fn_skf_z(1))
            endif
            if(myrank_skf_z(1) == 0) write(id_skf_z(1),'(A9,2(A15))') 'nt', 'u_tau', 'cf'
        endif
        if (out_skf_z(2)) then
            call initCommForSkf(nz, st(3), en(3), comm_skf_z(2), myrank_skf_z(2))
            if (is_restart) then
                open(unit=id_skf_z(2), file=fn_skf_z(2), position='append')
            else
                open(unit=id_skf_z(2), file=fn_skf_z(2))
            endif
            if(myrank_skf_z(2) == 0) write(id_skf_z(2),'(A9,2(A15))') 'nt', 'u_tau', 'cf'
        endif

        if (out_region) then
            ! exclude the points out of the domain
            region_ijk(1) = max(region_ijk(1), 1)
            region_ijk(2) = min(region_ijk(2), nx)
            region_ijk(3) = max(region_ijk(3), 1)
            region_ijk(4) = min(region_ijk(4), ny)
            region_ijk(5) = max(region_ijk(5), 1)
            region_ijk(6) = min(region_ijk(6), nz)
            
            ext_region%have_this_region = .false.
            ext_region%gsize(1) = region_ijk(2) - region_ijk(1) + 1
            ext_region%gsize(2) = region_ijk(4) - region_ijk(3) + 1
            ext_region%gsize(3) = region_ijk(6) - region_ijk(5) + 1
            ext_region%stg = 1
            ext_region%eng = 1
            ext_region%offset = 0
            ext_region%st = 1
            ext_region%sz = 1

            if (region_ijk(1)<=en(1) .and. region_ijk(2)>=st(1) .and. &
                region_ijk(3)<=en(2) .and. region_ijk(4)>=st(2) .and. &
                region_ijk(5)<=en(3) .and. region_ijk(6)>=st(3)) then
                
                ext_region%have_this_region = .true.
                ext_region%stg(:) = max(region_ijk(1:5:2), st(:))
                ext_region%eng(:) = min(region_ijk(2:6:2), en(:))
                ext_region%offset(:) = ext_region%stg(:) - region_ijk(1:5:2)
                ext_region%st(:) = ext_region%stg(:) - st(:) + 1
                ext_region%sz(:) = ext_region%eng(:) - ext_region%stg(:) + 1
                
                allocate( ext_buffer(ext_region%sz(1), ext_region%sz(2), ext_region%sz(3)) , stat=istatus)
                if (istatus /= 0) call abort(106, "initMonitor: Out of memory when allocating buffer for extracting region!")
                
            endif
        endif

        return
    end subroutine initMonitor

    subroutine outputMonitor(nt, u, v, w, p, u_crf)
        implicit none
        integer, intent(in) :: nt
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        real(fp), dimension(0:,0:,0:), intent(in) :: p
        real(fp), intent(in) :: u_crf

        if (out_forcing) call calcAndWriteForcing(id_forcing, nt, u, v, w, u_crf)
        
        if (out_probe_point) call extractAndWriteProbePoint(id_probe, nt, u, v, w, p, u_crf)

        if (out_skf_z(1)) call calcAndWriteSkf(id_skf_z(1),comm_skf_z(1),myrank_skf_z(1),nt, &
                                               1,0.5*dzf( 1),u, u_crf)
        if (out_skf_z(2)) call calcAndWriteSkf(id_skf_z(2),comm_skf_z(2),myrank_skf_z(2),nt, &
                                               sz(3),0.5*dzf(sz(3)),u, u_crf)
        if ( out_region .and. mod(nt, nt_out_region)==0 ) call extractAndWriteRegion(nt, u, v, w, p, u_crf)

        return
    end subroutine outputMonitor

    subroutine freeMonitor()
        implicit none

        if (out_forcing) close(id_forcing)
        if (out_probe_point) close(id_probe)
        if (out_skf_z(1)) close(id_skf_z(1))
        if (out_skf_z(2)) close(id_skf_z(2))
        if (out_region .and. allocated(ext_buffer)) deallocate(ext_buffer)
        
        return
    end subroutine freeMonitor

    subroutine extractAndWriteProbePoint(id_out, nt, u, v, w, p, u_crf)
        implicit none
        integer, intent(in) :: id_out
        integer, intent(in) :: nt
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        real(fp), dimension(0:,0:,0:), intent(in) :: p
        real(fp), intent(in) :: u_crf

        integer :: i, j, k

        if (probe_point%have_this_point) then
            i = probe_point%i
            j = probe_point%j
            k = probe_point%k
            probe_point%u = u(i, j, k) + u_crf
            probe_point%v = v(i, j, k)
            probe_point%w = w(i, j, k)
            probe_point%p = p(i, j, k)

            write(id_out,'(I9,4(E15.7))') nt, probe_point%u, probe_point%v, probe_point%w, probe_point%p 
        endif

        return
    end subroutine extractAndWriteProbePoint

    subroutine initCommForSkf(k_wall, ks, ke, comm_new, myrank_new)
        implicit none
        integer, intent(in) :: k_wall, ks, ke
        integer, intent(out) :: comm_new, myrank_new

        integer :: key, color

        key = 0
        color = MPI_UNDEFINED
        if (k_wall>=ks .and. k_wall<=ke) color = 1
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, key, comm_new, ierr)
        if (comm_new == MPI_COMM_NULL) then
            myrank_new = myrank
        else
            call MPI_COMM_RANK(comm_new, myrank_new, ierr)
        endif

        return
    end subroutine initCommForSkf

    subroutine calcAndWriteSkf(id_out, comm_skf, myrank_skf, nt, k_wall, dz, vel, vel_crf)
        implicit none
        integer, intent(in) :: id_out
        integer, intent(in) :: comm_skf, myrank_skf
        integer, intent(in) :: nt, k_wall
        real(fp), intent(in) :: dz
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: vel
        real(fp), intent(in) :: vel_crf

        real(fp) :: my_tau_w, tau_w, u_tau, cf
        integer :: i, j

        if (comm_skf /= MPI_COMM_NULL) then
            my_tau_w = 0.0_fp
            !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) REDUCTION(+:my_tau_w) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j)
            do j = 1, sz(2)
            do i = 1, sz(1)
                my_tau_w = my_tau_w + re_inv*( vel(i,j,k_wall)+vel_crf )/dz
            enddo
            enddo
            !$OMP END PARALLEL DO
            call MPI_REDUCE(my_tau_w, tau_w, 1, MPI_REAL_FP, MPI_SUM, 0, comm_skf, ierr)
            if (myrank_skf == 0) then
                tau_w = tau_w/nx/ny
                u_tau = sqrt(tau_w)
                cf = tau_w/(0.5_fp*u_ref*u_ref)

                write(id_out,'(I9,2(E15.7))') nt, u_tau, cf
            endif
        endif
        
        return
    end subroutine calcAndWriteSkf

    subroutine calcAndWriteForcing(id_out, nt, u, v, w, u_crf)
        implicit none
        integer, intent(in) :: id_out
        integer, intent(in) :: nt
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        real(fp), intent(in) :: u_crf

        real(fp), dimension(3) :: vel_mean
        
        if (any(is_forced(:))) then
            vel_mean(:) = 0.0_fp
            if (is_forced(1)) vel_mean(1) = Mean_h(nhalo, sz, nx, ny, dzflzi, u)
            if (is_forced(2)) vel_mean(2) = Mean_h(nhalo, sz, nx, ny, dzflzi, v)
            if (is_forced(3)) vel_mean(3) = Mean_h(nhalo, sz, nx, ny, dzflzi, w)
            
            if (myrank == 0) then
                write(id_out,'(I9,3(E15.7))') nt, vel_mean(1)+u_crf, vel_mean(2), vel_mean(3)
            endif
        endif

        return
    end subroutine calcAndWriteForcing

    subroutine extractAndWriteRegion(nt, u, v, w, p, u_crf)
        implicit none
        integer, intent(in) :: nt
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        real(fp), dimension(0:,0:,0:), intent(in) :: p
        real(fp), intent(in) :: u_crf

        character(12) :: string_dump
        integer(HID_T) :: fh

        write(string_dump, "('_',I8.8,'.h5')") nt
        call createFile(trim(adjustl(fn_prefix_region))//string_dump, fh)
        call writeAttribute(fh, 'nt', nt)

        if (ext_region%have_this_region) then
            call extractRegion(nhalo, u, 1, ext_region%st, ext_region%sz, ext_buffer)
            ext_buffer(:,:,:) = ext_buffer(:,:,:) + u_crf
        endif
        call write3d(fh, ext_region%have_this_region, ext_region%gsize, ext_region%offset+1, ext_region%sz, &
                     'u', ext_buffer)
        if (ext_region%have_this_region) call extractRegion(nhalo, v, 2, ext_region%st, ext_region%sz, ext_buffer)
        call write3d(fh, ext_region%have_this_region, ext_region%gsize, ext_region%offset+1, ext_region%sz, &
                     'v', ext_buffer)
        if (ext_region%have_this_region) call extractRegion(nhalo, w, 3, ext_region%st, ext_region%sz, ext_buffer)
        call write3d(fh, ext_region%have_this_region, ext_region%gsize, ext_region%offset+1, ext_region%sz, &
                     'w', ext_buffer)
        if (ext_region%have_this_region) call extractRegion((/1,1,1,1,1,1/), p, 0, ext_region%st, ext_region%sz, ext_buffer)
        call write3d(fh, ext_region%have_this_region, ext_region%gsize, ext_region%offset+1, ext_region%sz, &
                     'p', ext_buffer)
        
        call closeFile(fh)

        return
    end subroutine extractAndWriteRegion

    subroutine extractRegion(nhalo, var, offset_dir, reg_st, reg_sz, buffer)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: var
        integer, intent(in) :: offset_dir
        integer, dimension(3), intent(in) :: reg_st, reg_sz
        real(fp), dimension(:,:,:), intent(out) :: buffer

        integer :: i, j, k
        integer :: is, js, ks

        is = reg_st(1)-1
        js = reg_st(2)-1
        ks = reg_st(3)-1

        select case(offset_dir)
        case(0)
            !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
            do k = 1, reg_sz(3)
            do j = 1, reg_sz(2)
            do i = 1, reg_sz(1)
                buffer(i, j, k) = var(i+is, j+js, k+ks)
            enddo
            enddo
            enddo
        case(1)
            !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
            do k = 1, reg_sz(3)
            do j = 1, reg_sz(2)
            do i = 1, reg_sz(1)
                buffer(i, j, k) = 0.5*( var(i-1+is, j+js, k+ks) + var(i+is, j+js, k+ks) )
            enddo
            enddo
            enddo
        case(2)
            !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
            do k = 1, reg_sz(3)
            do j = 1, reg_sz(2)
            do i = 1, reg_sz(1)
                buffer(i, j, k) = 0.5*( var(i+is, j-1+js, k+ks) + var(i+is, j+js, k+ks) )
            enddo
            enddo
            enddo
        case(3)
            !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
            do k = 1, reg_sz(3)
            do j = 1, reg_sz(2)
            do i = 1, reg_sz(1)
                buffer(i, j, k) = 0.5*( var(i+is, j+js, k-1+ks) + var(i+is, j+js, k+ks) )
            enddo
            enddo
            enddo
        end select

        return
    end subroutine extractRegion

end module mod_monitor
