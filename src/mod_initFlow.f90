module mod_initFlow
    use mod_type,       only: fp
    use mod_parameters, only: re, initial_field, u_ref, u0, init_with_noise, noise_intensity, &
                              lz, nx, ny, nz, nhalo
    use mod_mpi,        only: sz, st, en, MPI_REAL_FP, ierr, myrank
    use mod_mesh,       only: zc_global, dzflzi
    use mod_utils,      only: Mean_h
    !$ use omp_lib
    
    implicit none
    ! make everything private unless declared public
    private

    ! public user routines
    public :: initFlow
    
contains

    subroutine initFlow(u, v, w)
        implicit none
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: u, v, w

        real(fp), dimension(1:sz(3)) :: zclzi, u_prof
        real(fp) :: meanold
        logical :: is_mean
        integer :: i, j, k

        select case(trim(adjustl(initial_field)))
        case('uni')
            u_prof(:) = u0
            is_mean = .false.
        case('poi')
            zclzi(:) = zc_global(st(3):en(3))/lz
            call setProfile_Poiseuille(sz(3), zclzi, u_ref, u_prof)
            is_mean = .true.
        case('log')
            zclzi(:) = zc_global(st(3):en(3))/lz
            call setProfile_Log(sz(3), zclzi, re, u_prof)
            is_mean = .true.
        case('mpi')
            u_prof(:) = real(myrank,fp)
        end select

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            u (i,j,k) = u_prof(k);
            ! v (i,j,k) = 0.0_fp     ! default value
            ! w (i,j,k) = 0.0_fp     ! default value
            ! p (i,j,k) = 0.0_fp     ! default value
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        if (init_with_noise) then
            call addNoise(st,sz,(/nx,ny,nz/),123,noise_intensity,u(1:sz(1),1:sz(2),1:sz(3)))
            call addNoise(st,sz,(/nx,ny,nz/),456,noise_intensity,v(1:sz(1),1:sz(2),1:sz(3)))
            call addNoise(st,sz,(/nx,ny,nz/),789,noise_intensity,w(1:sz(1),1:sz(2),1:sz(3)))
        endif

        if (is_mean) then
            meanold = Mean_h(nhalo, sz, nx, ny, dzflzi, u)
            if(meanold.ne.0.0_fp) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) &
                !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
                do k = 1, sz(3)
                do j = 1, sz(2)
                do i = 1, sz(1)
                    u(i,j,k) = u(i,j,k)/meanold*u_ref
                enddo
                enddo
                enddo
                !$OMP END PARALLEL DO
            endif
        endif

        return
    end subroutine initFlow

    subroutine setProfile_Poiseuille(n, z_norm, factor, vel)
        implicit none
        integer, intent(in) :: n
        real(fp), dimension(:), intent(in) :: z_norm
        real(fp), intent(in) :: factor
        real(fp), dimension(:), intent(out) :: vel

        integer :: k

        do k = 1, n
            vel(k) = 6.0_fp*(1.0_fp-z_norm(k))*z_norm(k) * factor
        enddo

        return
    end subroutine setProfile_Poiseuille

    subroutine setProfile_Log(n, z_norm, re, vel)
        implicit none
        integer, intent(in) :: n
        real(fp), dimension(:), intent(in) :: z_norm
        real(fp), intent(in) :: re
        real(fp), dimension(:), intent(out) :: vel

        real(fp) :: retau, z
        integer :: k

        retau = 0.09_fp*re**(0.88_fp)
        do k = 1, n
            z = min(z_norm(k),1.0_fp-z_norm(k))*2.0_fp*retau
            if (z <= 11.6) then
                vel(k) = z
            else
                vel(k) = 2.5_fp*log(z) + 5.5_fp
            endif
        enddo

        return
    end subroutine setProfile_Log

    subroutine addNoise(st,sz,sz_global,iseed,norm,var)
        implicit none
        integer, dimension(3), intent(in) :: st, sz, sz_global
        integer , intent(in) :: iseed
        real(fp), intent(in) :: norm 
        real(fp), dimension(:,:,:), intent(inout) :: var

        integer(4), dimension(64) :: seed
        real(fp) :: rn
        integer :: i,j,k,ii,jj,kk

        ! seed(:) = iseed
        ! call random_seed( put = seed )
        ! do k = 1, sz_global(3)
        !     kk = k-st(3)+1
        !     do j = 1, sz_global(2)
        !         jj = j-st(2)+1
        !         do i = 1, sz_global(1)
        !             ii = i-st(1)+1
        !             call random_number(rn)
        !             if(ii>=1.and.ii<=sz(1) .and. &
        !                jj>=1.and.jj<=sz(2) .and. &
        !                kk>=1.and.kk<=sz(3)) then
        !                 var(ii,jj,kk) = var(ii,jj,kk) + 2.0_fp*(rn-0.5_fp)*norm
        !             endif
        !         enddo
        !     enddo
        ! enddo

        seed(:) = iseed + myrank
        call random_seed( put = seed )
        do k = 1, sz(3)
            do j = 1, sz(2)
                do i = 1, sz(1)
                    call random_number(rn)
                    var(i,j,k) = var(i,j,k) + 2.0_fp*(rn-0.5_fp)*norm
                enddo
            enddo
        enddo

        return
    end subroutine addNoise

end module mod_initFlow