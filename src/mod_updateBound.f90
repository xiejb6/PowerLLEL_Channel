module mod_updateBound
    use mod_type,       only: fp
    use mod_parameters, only: nhalo, nhalo_one
    use mod_mpi,        only: sz, halotype_vel, halotype_one, neighbor, comm_cart, ierr
#ifdef NB_HALO
    use mod_mpi,        only: mpitype_nbhalo_vel, neighbor_nbhalo
#endif
    !$ use omp_lib
#ifdef GPTL
    use gptl
#endif

    implicit none

    include 'mpif.h'

    ! make everything private unless declared public
    private

#ifdef GPTL
    integer :: ret
#endif

    public :: updateBoundVel, imposeBCVel, updateBoundCenteredVel, updateBoundP
#ifdef NB_HALO
    public :: updateHaloISend, updateHaloIRecv, updateHaloWaitall
#endif
    
contains
    subroutine updateBoundVel(u, v, w, u_crf)
        implicit none
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: u, v, w
        real(fp), intent(in) :: u_crf

        integer :: ibound, idir

#ifdef GPTL
        ret = gptlstart('--Update halo vel')
#endif

        call updateHalo(nhalo, halotype_vel, u)
        call updateHalo(nhalo, halotype_vel, v)
        call updateHalo(nhalo, halotype_vel, w)

#ifdef GPTL
        ret = gptlstop('--Update halo vel')
        ret = gptlstart('--Impose BC vel')
#endif
        call imposeBCVel(u, v, w, u_crf)
#ifdef GPTL
        ret = gptlstop('--Impose BC vel')
#endif
        
        return
    end subroutine updateBoundVel

    subroutine imposeBCVel(u, v, w, u_crf)
        implicit none
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: u, v, w
        real(fp), intent(in) :: u_crf

        integer :: ibound, idir
        
        ! B.C. in x direction
        idir = 1
        if (neighbor(1) == MPI_PROC_NULL) then
            ibound = 0
            call imposePeriodicBC(nhalo, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo, sz, ibound, idir, w)
        endif
        if (neighbor(2) == MPI_PROC_NULL) then
            ibound = 1
            call imposePeriodicBC(nhalo, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo, sz, ibound, idir, w)
        endif
        ! B.C. in y direction
        idir = 2
        if (neighbor(3) == MPI_PROC_NULL) then
            ibound = 0
            call imposePeriodicBC(nhalo, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo, sz, ibound, idir, w)
        endif
        if (neighbor(4) == MPI_PROC_NULL) then
            ibound = 1
            call imposePeriodicBC(nhalo, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo, sz, ibound, idir, w)
        endif
        ! B.C. in z direction
        idir = 3
        if (neighbor(5) == MPI_PROC_NULL) then
            ibound = 0
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  u, u_crf)
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  v)
            call imposeNoSlipBC(nhalo, sz, ibound, .false., w)
        endif
        if (neighbor(6) == MPI_PROC_NULL) then
            ibound = 1
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  u, u_crf)
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  v)
            call imposeNoSlipBC(nhalo, sz, ibound, .false., w)
        endif

        return
    end subroutine imposeBCVel

    subroutine updateBoundCenteredVel(u, v, w)
        implicit none
        real(fp), dimension(0:,0:,0:), intent(inout) :: u, v, w

        integer :: ibound, idir
        integer :: i, j, k

        call updateHalo(nhalo_one, halotype_one, u)
        call updateHalo(nhalo_one, halotype_one, v)
        call updateHalo(nhalo_one, halotype_one, w)

        ! B.C. in x direction
        idir = 1
        if (neighbor(1) == MPI_PROC_NULL) then
            ibound = 0
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, w)
        endif
        if (neighbor(2) == MPI_PROC_NULL) then
            ibound = 1
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, w)
        endif
        ! B.C. in y direction
        idir = 2
        if (neighbor(3) == MPI_PROC_NULL) then
            ibound = 0
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, w)
        endif
        if (neighbor(4) == MPI_PROC_NULL) then
            ibound = 1
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, u)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, v)
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, w)
        endif
        ! B.C. in z direction
        idir = 3
        if (neighbor(5) == MPI_PROC_NULL) then
            !$OMP PARALLEL DO SCHEDULE(STATIC)
            do j = 1-nhalo_one(3), sz(2)+nhalo_one(4)
            do i = 1-nhalo_one(1), sz(1)+nhalo_one(2)
                u(i,j,0) = 0.0_fp
                v(i,j,0) = 0.0_fp
                w(i,j,0) = 0.0_fp
            enddo
            enddo
            !$OMP END PARALLEL DO
            ! ibound = 0
            ! call imposeNoSlipBC(nhalo_one, sz, ibound, .false., u)
            ! call imposeNoSlipBC(nhalo_one, sz, ibound, .false., v)
            ! call imposeNoSlipBC(nhalo_one, sz, ibound, .false., w)
        endif
        if (neighbor(6) == MPI_PROC_NULL) then
            k = sz(3)+1
            !$OMP PARALLEL DO SCHEDULE(STATIC)
            do j = 1-nhalo_one(3), sz(2)+nhalo_one(4)
            do i = 1-nhalo_one(1), sz(1)+nhalo_one(2)
                u(i,j,k) = 0.0_fp
                v(i,j,k) = 0.0_fp
                w(i,j,k) = 0.0_fp
            enddo
            enddo
            !$OMP END PARALLEL DO
            ! ibound = 1
            ! call imposeNoSlipBC(nhalo_one, sz, ibound, .true., u)
            ! call imposeNoSlipBC(nhalo_one, sz, ibound, .true., v)
            ! call imposeNoSlipBC(nhalo_one, sz, ibound, .true., w)
        endif
        
        return
    end subroutine updateBoundCenteredVel
    
    subroutine updateBoundP(p)
        implicit none
        real(fp), dimension(0:,0:,0:), intent(inout) :: p

        integer :: ibound, idir

#ifdef GPTL
        ret = gptlstart('--Update halo pres')
#endif

        call updateHalo(nhalo_one, halotype_one, p)

#ifdef GPTL
        ret = gptlstop('--Update halo pres')
        ret = gptlstart('--Impose BC pres')
#endif

        ! B.C. in x direction
        idir = 1
        if (neighbor(1) == MPI_PROC_NULL) then
            ibound = 0
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, p)
        endif
        if (neighbor(2) == MPI_PROC_NULL) then
            ibound = 1
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, p)
        endif
        ! B.C. in y direction
        idir = 2
        if (neighbor(3) == MPI_PROC_NULL) then
            ibound = 0
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, p)
        endif
        if (neighbor(4) == MPI_PROC_NULL) then
            ibound = 1
            call imposePeriodicBC(nhalo_one, sz, ibound, idir, p)
        endif
        ! B.C. in z direction
        idir = 3
        if (neighbor(5) == MPI_PROC_NULL) then
            ibound = 0
            call imposeZeroGradBC(nhalo_one, sz, ibound, p)
        endif
        if (neighbor(6) == MPI_PROC_NULL) then
            ibound = 1
            call imposeZeroGradBC(nhalo_one, sz, ibound, p)
        endif

#ifdef GPTL
        ret = gptlstop('--Impose BC pres')
#endif
        
        return
    end subroutine updateBoundP

    subroutine imposePeriodicBC(nhalo, sz, ibound, idir, var)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        integer, intent(in) :: ibound, idir
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var

        integer :: i, j, k

        if (idir == 1) then
            select case(ibound)
            case(0)
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do k = 1-nhalo(5), sz(3)+nhalo(6)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                    var(1-nhalo(1):0, j, k) = var(sz(1)+1-nhalo(1):sz(1), j, k)
                enddo
                enddo
                !$OMP END PARALLEL DO
            case(1)
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do k = 1-nhalo(5), sz(3)+nhalo(6)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                    var(sz(1)+1:sz(1)+nhalo(2), j, k) = var(1:nhalo(2), j, k)
                enddo
                enddo
                !$OMP END PARALLEL DO
            end select
        endif

        return
    end subroutine imposePeriodicBC

    subroutine imposeNoSlipBC(nhalo, sz, ibound, centered, var, vel_crf)
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        integer, intent(in) :: ibound
        logical, intent(in) :: centered
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var
        real(fp), optional, intent(in) :: vel_crf

        real(fp) :: bcvalue
        integer :: i, j, k, n

        bcvalue = 0.0_fp
        if (present(vel_crf)) bcvalue = 0.0_fp - vel_crf

        select case(ibound)
        case(0)
            if (centered) then
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,0  ) = 2.0_fp*bcvalue - var(i,j,1)
                enddo
                enddo
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,0  ) = bcvalue
                enddo
                enddo
                !$OMP END PARALLEL DO
            endif
        case(1)
            n = sz(3)
            if (centered) then
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,n+1) = 2.0_fp*bcvalue - var(i,j,n)
                enddo
                enddo
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,n  ) = bcvalue
                    var(i,j,n+1) = bcvalue
                enddo
                enddo
                !$OMP END PARALLEL DO
            endif
        end select

        return
    end subroutine imposeNoSlipBC

    subroutine imposeZeroGradBC(nhalo, sz, ibound, var)
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        integer, intent(in) :: ibound
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var

        integer :: i, j, k, n

        select case(ibound)
        case(0)
            !$OMP PARALLEL DO SCHEDULE(STATIC)
            do j = 1-nhalo(3), sz(2)+nhalo(4)
            do i = 1-nhalo(1), sz(1)+nhalo(2)
                var(i,j,0  ) = var(i,j,1)
            enddo
            enddo
            !$OMP END PARALLEL DO
        case(1)
            n = sz(3)
            !$OMP PARALLEL DO SCHEDULE(STATIC)
            do j = 1-nhalo(3), sz(2)+nhalo(4)
            do i = 1-nhalo(1), sz(1)+nhalo(2)
                var(i,j,n+1) = var(i,j,n)
            enddo
            enddo
            !$OMP END PARALLEL DO

        end select

        return
    end subroutine imposeZeroGradBC
    
    subroutine updateHalo(nhalo, halotype, var)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(6), intent(in) :: halotype
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var

        integer :: tag_halo = 0

        ! update halo cells
        ! *** west/east ***
        ! call MPI_SENDRECV(var(1      , 1-nhalo(3), 1-nhalo(5)), 1, halotype(2), neighbor(1), tag_halo, &
        !                   var(sz(1)+1, 1-nhalo(3), 1-nhalo(5)), 1, halotype(2), neighbor(2), tag_halo, &
        !                   comm_cart, MPI_STATUS_IGNORE, ierr)
        ! call MPI_SENDRECV(var(sz(1)+1-nhalo(1), 1-nhalo(3), 1-nhalo(5)), 1, halotype(1), neighbor(2), tag_halo, &
        !                   var(1-nhalo(1)      , 1-nhalo(3), 1-nhalo(5)), 1, halotype(1), neighbor(1), tag_halo, &
        !                   comm_cart, MPI_STATUS_IGNORE, ierr)
        ! *** south/north ***
        call MPI_SENDRECV(var(1-nhalo(1), 1      , 1-nhalo(5)), 1, halotype(4), neighbor(3), tag_halo, &
                          var(1-nhalo(1), sz(2)+1, 1-nhalo(5)), 1, halotype(4), neighbor(4), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        call MPI_SENDRECV(var(1-nhalo(1), sz(2)+1-nhalo(3), 1-nhalo(5)), 1, halotype(3), neighbor(4), tag_halo, &
                          var(1-nhalo(1), 1-nhalo(3)      , 1-nhalo(5)), 1, halotype(3), neighbor(3), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        ! *** bottom/top ***
        call MPI_SENDRECV(var(1-nhalo(1), 1-nhalo(3), 1      ), 1, halotype(6), neighbor(5), tag_halo, &
                          var(1-nhalo(1), 1-nhalo(3), sz(3)+1), 1, halotype(6), neighbor(6), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        call MPI_SENDRECV(var(1-nhalo(1), 1-nhalo(3), sz(3)+1-nhalo(5)), 1, halotype(5), neighbor(6), tag_halo, &
                          var(1-nhalo(1), 1-nhalo(3), 1-nhalo(5)      ), 1, halotype(5), neighbor(5), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)

        return
    end subroutine updateHalo

#ifdef NB_HALO

    subroutine updateHaloIRecv(nhalo, tag, var, irecv_req)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, intent(in) :: tag
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var
        integer, dimension(8), intent(out) :: irecv_req

        integer :: ist, jst, kst

        ist = 1-nhalo(1)
        jst = 1-nhalo(3)
        kst = 1-nhalo(5)

        call MPI_IRECV(var(ist, sz(2)+1, 1), 1, mpitype_nbhalo_vel(2), neighbor_nbhalo(2), &
                       tag, comm_cart, irecv_req(2), ierr)   ! y1 recv y0 send
        call MPI_IRECV(var(ist, jst    , 1), 1, mpitype_nbhalo_vel(1), neighbor_nbhalo(1), &
                       tag, comm_cart, irecv_req(1), ierr)   ! y0 recv y1 send
        call MPI_IRECV(var(ist, 1, sz(3)+1), 1, mpitype_nbhalo_vel(4), neighbor_nbhalo(4), &
                       tag, comm_cart, irecv_req(4), ierr)   ! z1 recv z0 send
        call MPI_IRECV(var(ist, 1, kst    ), 1, mpitype_nbhalo_vel(3), neighbor_nbhalo(3), &
                       tag, comm_cart, irecv_req(3), ierr)   ! z0 recv z1 send
        call MPI_IRECV(var(ist, sz(2)+1, sz(3)+1), 1, mpitype_nbhalo_vel(8), neighbor_nbhalo(8), &
                       tag, comm_cart, irecv_req(8), ierr) ! north_top recv south_bottom send
        call MPI_IRECV(var(ist,     jst,     kst), 1, mpitype_nbhalo_vel(5), neighbor_nbhalo(5), &
                       tag, comm_cart, irecv_req(5), ierr) ! south_bottom recv north_top send
        call MPI_IRECV(var(ist, sz(2)+1,     kst), 1, mpitype_nbhalo_vel(7), neighbor_nbhalo(7), &
                       tag, comm_cart, irecv_req(7), ierr) ! north_bottom recv south_top send
        call MPI_IRECV(var(ist,     jst, sz(3)+1), 1, mpitype_nbhalo_vel(6), neighbor_nbhalo(6), &
                       tag, comm_cart, irecv_req(6), ierr) ! south_top recv north_bottom send

        return
    end subroutine updateHaloIRecv

    subroutine updateHaloISend(nhalo, tag, var, isend_req)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, intent(in) :: tag
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var
        integer, dimension(8), intent(out) :: isend_req

        integer :: ist, jst, kst

        ist = 1-nhalo(1)
        jst = 1-nhalo(3)
        kst = 1-nhalo(5)

        call MPI_ISEND(var(ist,         1, 1), 1, mpitype_nbhalo_vel(2), neighbor_nbhalo(1), &
                       tag, comm_cart, isend_req(1), ierr)   ! y1 recv y0 send
        call MPI_ISEND(var(ist, sz(2)+jst, 1), 1, mpitype_nbhalo_vel(1), neighbor_nbhalo(2), &
                       tag, comm_cart, isend_req(2), ierr)   ! y0 recv y1 send
        call MPI_ISEND(var(ist, 1,         1), 1, mpitype_nbhalo_vel(4), neighbor_nbhalo(3), &
                       tag, comm_cart, isend_req(3), ierr)   ! z1 recv z0 send
        call MPI_ISEND(var(ist, 1, sz(3)+kst), 1, mpitype_nbhalo_vel(3), neighbor_nbhalo(4), &
                       tag, comm_cart, isend_req(4), ierr)   ! z0 recv z1 send
        call MPI_ISEND(var(ist,         1,         1), 1, mpitype_nbhalo_vel(8), neighbor_nbhalo(5), &
                       tag, comm_cart, isend_req(5), ierr) ! north_top recv south_bottom send
        call MPI_ISEND(var(ist, sz(2)+jst, sz(3)+kst), 1, mpitype_nbhalo_vel(5), neighbor_nbhalo(8), &
                       tag, comm_cart, isend_req(8), ierr) ! south_bottom recv north_top send
        call MPI_ISEND(var(ist,         1, sz(3)+kst), 1, mpitype_nbhalo_vel(7), neighbor_nbhalo(6), &
                       tag, comm_cart, isend_req(6), ierr) ! north_bottom recv south_top send
        call MPI_ISEND(var(ist, sz(2)+jst,         1), 1, mpitype_nbhalo_vel(6), neighbor_nbhalo(7), &
                       tag, comm_cart, isend_req(7), ierr) ! south_top recv north_bottom send

        return
    end subroutine updateHaloISend

    subroutine updateHaloWaitall(isend_req, irecv_req)
        implicit none
        integer, dimension(:), intent(inout) :: isend_req, irecv_req

        call MPI_WAITALL(size(isend_req), isend_req, MPI_STATUS_IGNORE, ierr)
        call MPI_WAITALL(size(irecv_req), irecv_req, MPI_STATUS_IGNORE, ierr)

        return
    end subroutine updateHaloWaitall
    
    subroutine updateNBHalo(nhalo, halotype, var)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(8), intent(in) :: halotype
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var

        integer :: tag_halo = 0

        ! update halo cells
        ! *** south/north ***
        call MPI_SENDRECV(var(1-nhalo(1),       1, 1), 1, halotype(2), neighbor_nbhalo(1), tag_halo, &
                          var(1-nhalo(1), sz(2)+1, 1), 1, halotype(2), neighbor_nbhalo(2), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        call MPI_SENDRECV(var(1-nhalo(1), sz(2)+1-nhalo(3), 1), 1, halotype(1), neighbor_nbhalo(2), tag_halo, &
                          var(1-nhalo(1),       1-nhalo(3), 1), 1, halotype(1), neighbor_nbhalo(1), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        ! *** bottom/top ***
        call MPI_SENDRECV(var(1-nhalo(1), 1,       1), 1, halotype(4), neighbor_nbhalo(3), tag_halo, &
                          var(1-nhalo(1), 1, sz(3)+1), 1, halotype(4), neighbor_nbhalo(4), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        call MPI_SENDRECV(var(1-nhalo(1), 1, sz(3)+1-nhalo(5)), 1, halotype(3), neighbor_nbhalo(4), tag_halo, &
                          var(1-nhalo(1), 1,       1-nhalo(5)), 1, halotype(3), neighbor_nbhalo(3), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        ! *** south_bottom/north_top ***
        call MPI_SENDRECV(var(1-nhalo(1),       1,       1), 1, halotype(8), neighbor_nbhalo(5), tag_halo, &
                          var(1-nhalo(1), sz(2)+1, sz(3)+1), 1, halotype(8), neighbor_nbhalo(8), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        call MPI_SENDRECV(var(1-nhalo(1), sz(2)+1-nhalo(3), sz(3)+1-nhalo(5)), 1, halotype(5), neighbor_nbhalo(8), tag_halo, &
                        var(1-nhalo(1),       1-nhalo(3),       1-nhalo(5)), 1, halotype(5), neighbor_nbhalo(5), tag_halo, &
                        comm_cart, MPI_STATUS_IGNORE, ierr)
        ! *** south_top/north_bottom ***
        call MPI_SENDRECV(var(1-nhalo(1),       1, sz(3)+1-nhalo(5)), 1, halotype(7), neighbor_nbhalo(6), tag_halo, &
                          var(1-nhalo(1), sz(2)+1,       1-nhalo(5)), 1, halotype(7), neighbor_nbhalo(7), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        call MPI_SENDRECV(var(1-nhalo(1), sz(2)+1-nhalo(3),       1), 1, halotype(6), neighbor_nbhalo(7), tag_halo, &
                          var(1-nhalo(1),       1-nhalo(3), sz(3)+1), 1, halotype(6), neighbor_nbhalo(6), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)

        return
    end subroutine updateNBHalo

#endif

end module mod_updateBound