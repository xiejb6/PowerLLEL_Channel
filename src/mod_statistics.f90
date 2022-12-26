module mod_statistics
    use mod_type
    use mod_parameters, only: nt_init_stat, stat_info, stat_which_var, sample_interval, nhalo
    use mod_utils,      only: abort, setZero
    !$ use omp_lib
    
    implicit none
    ! make everything private unless declared public
    private

    real(fp), allocatable, dimension(:,:,:), save, public :: u_stat,  v_stat,  w_stat
    real(fp), allocatable, dimension(:,:,:), save, public :: u2_stat, v2_stat, w2_stat
    real(fp), allocatable, dimension(:,:,:), save, public :: uv_stat, uw_stat, vw_stat
    real(fp), allocatable, dimension(:,:,:), save, public :: p_stat,  p2_stat

    integer, save :: ie, je, ke

    public :: allocStat, freeStat, initStat, calcStat

contains
    subroutine allocStat(sz)
        implicit none
        integer, dimension(3), intent(in) :: sz

        integer :: istatus

        ie = sz(1)
        je = sz(2)
        ke = sz(3)

        if (stat_which_var( 1)) then
            allocate(u_stat (1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 2)) then
            allocate(v_stat (1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 3)) then
            allocate(w_stat (1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 4)) then
            allocate(u2_stat(1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 5)) then
            allocate(v2_stat(1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 6)) then
            allocate(w2_stat(1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 7)) then
            allocate(uv_stat(1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 8)) then
            allocate(uw_stat(1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var( 9)) then
            allocate(vw_stat(1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var(10)) then
            allocate(p_stat (1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif
        if (stat_which_var(11)) then
            allocate(p2_stat(1:ie, 1:je, 1:ke), stat=istatus)
            if (istatus /= 0) call abort(102, "allocStat: Out of memory when allocating statistical variables!")
        endif

        return
    end subroutine allocStat

    subroutine freeStat()
        implicit none

        if (allocated(u_stat )) deallocate(u_stat )
        if (allocated(v_stat )) deallocate(v_stat )
        if (allocated(w_stat )) deallocate(w_stat )
        if (allocated(u2_stat)) deallocate(u2_stat)
        if (allocated(v2_stat)) deallocate(v2_stat)
        if (allocated(w2_stat)) deallocate(w2_stat)
        if (allocated(uv_stat)) deallocate(uv_stat)
        if (allocated(uw_stat)) deallocate(uw_stat)
        if (allocated(vw_stat)) deallocate(vw_stat)
        if (allocated(p_stat )) deallocate(p_stat )
        if (allocated(p2_stat)) deallocate(p2_stat)

        return
    end subroutine freeStat

    subroutine initStat()
        implicit none

        stat_info%nts  = nt_init_stat+1
        stat_info%nte  = nt_init_stat
        stat_info%nspl = 0

        if (stat_which_var( 1)) call setZero(u_stat)
        if (stat_which_var( 2)) call setZero(v_stat)
        if (stat_which_var( 3)) call setZero(w_stat)
        if (stat_which_var( 4)) call setZero(u2_stat)
        if (stat_which_var( 5)) call setZero(v2_stat)
        if (stat_which_var( 6)) call setZero(w2_stat)
        if (stat_which_var( 7)) call setZero(uv_stat)
        if (stat_which_var( 8)) call setZero(uw_stat)
        if (stat_which_var( 9)) call setZero(vw_stat)
        if (stat_which_var(10)) call setZero(p_stat)
        if (stat_which_var(11)) call setZero(p2_stat)

        return        
    end subroutine initStat

    subroutine calcStat(nt, u, v, w, p, u_crf)
        implicit none
        integer,  intent(in) :: nt
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        real(fp), dimension(0:,0:,0:), intent(in) :: p
        real(fp), intent(in) :: u_crf

        real(fp), dimension(1:ie,1:je) :: uc, vc, wc
        integer :: i, j, k

        if (mod(nt-nt_init_stat, sample_interval) == 0) then

            stat_info%nte  = stat_info%nte  + sample_interval
            stat_info%nspl = stat_info%nspl + 1

            !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, uc, vc, wc)
            !$OMP DO SCHEDULE(STATIC)
            do k = 1, ke
                do j = 1, je
                do i = 1, ie
                    uc(i, j) = (u(i-1,j,k)+u(i,j,k))*0.5_fp + u_crf
                    vc(i, j) = (v(i,j-1,k)+v(i,j,k))*0.5_fp
                    wc(i, j) = (w(i,j,k-1)+w(i,j,k))*0.5_fp
                enddo
                enddo
                if (stat_which_var( 1)) then
                    u_stat (:,:,k) = u_stat (:,:,k) + uc(:,:)
                endif
                if (stat_which_var( 2)) then
                    v_stat (:,:,k) = v_stat (:,:,k) + vc(:,:)
                endif
                if (stat_which_var( 3)) then
                    w_stat (:,:,k) = w_stat (:,:,k) + wc(:,:)
                endif
                if (stat_which_var( 4)) then
                    u2_stat(:,:,k) = u2_stat(:,:,k) + uc(:,:)*uc(:,:)
                endif
                if (stat_which_var( 5)) then
                    v2_stat(:,:,k) = v2_stat(:,:,k) + vc(:,:)*vc(:,:)
                endif
                if (stat_which_var( 6)) then
                    w2_stat(:,:,k) = w2_stat(:,:,k) + wc(:,:)*wc(:,:)
                endif
                if (stat_which_var( 7)) then
                    uv_stat(:,:,k) = uv_stat(:,:,k) + uc(:,:)*vc(:,:)
                endif
                if (stat_which_var( 8)) then
                    uw_stat(:,:,k) = uw_stat(:,:,k) + uc(:,:)*wc(:,:)
                endif                
                if (stat_which_var( 9)) then
                    vw_stat(:,:,k) = vw_stat(:,:,k) + vc(:,:)*wc(:,:)
                endif                
            enddo
            !$OMP END DO
            if (stat_which_var(10)) then
                !$OMP DO SCHEDULE(STATIC)
                do k = 1, ke
                    p_stat (:,:,k) = p_stat (:,:,k) + p(:,:,k)
                enddo
                !$OMP END DO
            endif
            if (stat_which_var(11)) then
                !$OMP DO SCHEDULE(STATIC)
                do k = 1, ke
                    p2_stat(:,:,k) = p2_stat(:,:,k) + p(:,:,k)*p(:,:,k)
                enddo
                !$OMP END DO
            endif
            !$OMP END PARALLEL

        endif

        return
    end subroutine calcStat

end module mod_statistics