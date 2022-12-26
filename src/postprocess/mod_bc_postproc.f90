module mod_bc_postproc
    use mod_type, only: fp
    use mod_mpi,  only: sz, halotype_vel, neighbor, comm_cart
    !$ use omp_lib
    
    implicit none

    include 'mpif.h'
    
    ! make everything private unless declared public
    private
    public :: updateBoundVel
    
contains
    subroutine updateBoundVel(nhalo, u, v, w)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: u, v, w

        integer :: ibound, idir

        call updateHalo(nhalo, halotype_vel, u)
        call updateHalo(nhalo, halotype_vel, v)
        call updateHalo(nhalo, halotype_vel, w)

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
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  u)
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  v)
            call imposeNoSlipBC(nhalo, sz, ibound, .false., w)
        endif
        if (neighbor(6) == MPI_PROC_NULL) then
            ibound = 1
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  u)
            call imposeNoSlipBC(nhalo, sz, ibound, .true.,  v)
            call imposeNoSlipBC(nhalo, sz, ibound, .false., w)
        endif
        
        return
    end subroutine updateBoundVel

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

    subroutine imposeNoSlipBC(nhalo, sz, ibound, centered, var)
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        integer, intent(in) :: ibound
        logical, intent(in) :: centered
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: var

        integer :: i, j, k

        select case(ibound)
        case(0)
            if (centered) then
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,0  ) = - var(i,j,1)
                enddo
                enddo
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,0  ) = 0.0_fp
                enddo
                enddo
                !$OMP END PARALLEL DO
            endif
        case(1)
            if (centered) then
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,sz(3)+1) = - var(i,j,sz(3))
                enddo
                enddo
                !$OMP END PARALLEL DO
            else
                !$OMP PARALLEL DO SCHEDULE(STATIC)
                do j = 1-nhalo(3), sz(2)+nhalo(4)
                do i = 1-nhalo(1), sz(1)+nhalo(2)
                    var(i,j,sz(3)  ) = 0.0_fp
                    var(i,j,sz(3)+1) = var(i,j,sz(3)-1)
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

        integer :: i, j, k

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
            !$OMP PARALLEL DO SCHEDULE(STATIC)
            do j = 1-nhalo(3), sz(2)+nhalo(4)
            do i = 1-nhalo(1), sz(1)+nhalo(2)
                var(i,j,sz(3)+1) = var(i,j,sz(3))
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
        integer :: ierr

        ! update halo cells
        ! *** west/east ***
        call MPI_SENDRECV(var(1      , 1-nhalo(3), 1-nhalo(5)), 1, halotype(2), neighbor(1), tag_halo, &
                          var(sz(1)+1, 1-nhalo(3), 1-nhalo(5)), 1, halotype(2), neighbor(2), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
        call MPI_SENDRECV(var(sz(1)+1-nhalo(1), 1-nhalo(3), 1-nhalo(5)), 1, halotype(1), neighbor(2), tag_halo, &
                          var(1-nhalo(1)      , 1-nhalo(3), 1-nhalo(5)), 1, halotype(1), neighbor(1), tag_halo, &
                          comm_cart, MPI_STATUS_IGNORE, ierr)
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

end module mod_bc_postproc