module mod_utils
    use mod_type, only: fp
    use mod_mpi, only: myrank, ierr, MPI_REAL_FP, st
    !$ use omp_lib
    
    implicit none

    include 'mpif.h'
    
    ! make everything private unless declared public
    private

    type value_index_pair_t
        real(fp) :: value
        integer :: rank
        integer :: ig, jg, kg
    end type
    integer, parameter :: blockcount = 5
    integer, dimension(blockcount), save :: blocklengths, displacements, types
    integer, save :: ValueIndexPair_type, MPI_MAX_INDEX

    interface setZero
        module procedure setZero_1d
        module procedure setZero_2d
        module procedure setZero_3d
    end interface setZero

    public :: value_index_pair_t, abort, setZero, Mean_h, calcMaxCFL, calcMaxDiv, &
              initCheckCFLAndDiv, freeCheckCFLAndDiv, checkCFL, checkDiv, checkNaN

contains

    subroutine abort(errcode, errmsg)
        implicit none
        integer, intent(in) :: errcode
        character(*), intent(in) :: errmsg

        ! errcode = 
        ! 101 - Out of memory when allocating major variables
        ! 102 - Out of memory when allocating statistical variables
        ! 103 - Out of memory when allocating coefficients of tridiagonal systems
        ! 104 - Out of memory when allocating work arrays for PDD algorithm
        ! 105 - Out of memory when allocating work arrays for transpositions
        ! 106 - Out of memory when allocating buffer for extracting region
        
        if (myrank == 0) then
            write(*,*) "PowerLLEL.ERROR."//errmsg
        endif
        call MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)

        return
    end subroutine abort

    subroutine setZero_1d(array)
        implicit none
        real(fp), dimension(:), intent(out) :: array

        integer :: i

        do i = lbound(array,1), ubound(array,1)
            array(i) = 0.0_fp
        enddo

        return
    end subroutine setZero_1d

    subroutine setZero_2d(array)
        implicit none
        real(fp), dimension(:,:), intent(out) :: array

        integer :: i, j

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j)
        do j = lbound(array,2), ubound(array,2)
        do i = lbound(array,1), ubound(array,1)
            array(i, j) = 0.0_fp
        enddo
        enddo
        !$OMP END PARALLEL DO

        return
    end subroutine setZero_2d
    
    subroutine setZero_3d(array)
        implicit none
        real(fp), dimension(:,:,:), intent(out) :: array

        integer :: i, j, k

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = lbound(array,3), ubound(array,3)
        do j = lbound(array,2), ubound(array,2)
        do i = lbound(array,1), ubound(array,1)
            array(i, j, k) = 0.0_fp
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        return
    end subroutine setZero_3d

    function Mean_h(nhalo, sz, nx_global, ny_global, dzflzi, var)
        implicit none
        real(fp) :: Mean_h
        integer, dimension(6) :: nhalo
        integer, dimension(3) :: sz
        integer :: nx_global, ny_global
        real(fp), dimension(0:) :: dzflzi
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):) :: var

        integer  :: i, j, k
        
        Mean_h = 0.0_fp
        !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:Mean_h) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            Mean_h = Mean_h + var(i, j, k)*dzflzi(k)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO
        call MPI_ALLREDUCE(MPI_IN_PLACE, Mean_h, 1, MPI_REAL_FP, MPI_SUM, MPI_COMM_WORLD, ierr)
        Mean_h = Mean_h/nx_global/ny_global

        return
    end function Mean_h

    subroutine initCheckCFLAndDiv(cfl_max, div_max)
        implicit none
        type(value_index_pair_t), intent(out) :: cfl_max, div_max

        integer :: i, disp

        blocklengths = (/1,1,1,1,1/)
        types = (/MPI_REAL_FP,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER,MPI_INTEGER/)
        displacements(1) = 0
        do i = 2, blockcount
            call MPI_TYPE_EXTENT(types(i-1), disp, ierr)
            displacements(i) = displacements(i-1) + disp
        enddo
        call MPI_TYPE_STRUCT(blockcount,blocklengths,displacements,types,ValueIndexPair_type,ierr)
        call MPI_TYPE_COMMIT(ValueIndexPair_type,ierr)
        call MPI_OP_CREATE(findMPIMaxIndex,.true.,MPI_MAX_INDEX,ierr)

        cfl_max%value = 0.0_fp
        cfl_max%rank = myrank
        cfl_max%ig = 1
        cfl_max%jg = 1
        cfl_max%kg = 1

        div_max%value = 0.0_fp
        div_max%rank = myrank
        div_max%ig = 1
        div_max%jg = 1
        div_max%kg = 1

        return
    end subroutine initCheckCFLAndDiv

    subroutine findMPIMaxIndex(invec, inoutvec, len, dtype)
        implicit none
        type(value_index_pair_t), intent(in) :: invec(*)
        type(value_index_pair_t), intent(inout) :: inoutvec(*)
        integer, intent(in) :: len
        integer, intent(in) :: dtype

        integer :: i

        do i = 1, len
            if (invec(i)%value > inoutvec(i)%value) then
                inoutvec(i) = invec(i)
            endif
        enddo

        return
    end subroutine findMPIMaxIndex

    subroutine freeCheckCFLAndDiv()
        implicit none

        call MPI_OP_FREE(MPI_MAX_INDEX, ierr)
        call MPI_TYPE_FREE(ValueIndexPair_type, ierr)

        return
    end subroutine freeCheckCFLAndDiv


    subroutine calcMaxCFL(nhalo, sz, dt, dxi, dyi, dzfi, u, v, w, cfl_max_root)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        real(fp), intent(in) :: dt, dxi, dyi
        real(fp), dimension(0:), intent(in) :: dzfi
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        type(value_index_pair_t), intent(out) :: cfl_max_root
        
        type(value_index_pair_t) :: cfl_max
        real(fp) :: cfl_tmp, cfl_max_thread
        integer :: i, j, k
        integer :: i_tmp, j_tmp, k_tmp

        cfl_max%value = 0.0_fp
        cfl_max%rank = myrank
        cfl_max%ig = 1
        cfl_max%jg = 1
        cfl_max%kg = 1
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP PRIVATE(i, j, k, cfl_tmp, cfl_max_thread, i_tmp, j_tmp, k_tmp)
        cfl_max_thread = cfl_max%value
        i_tmp = cfl_max%ig
        j_tmp = cfl_max%jg
        k_tmp = cfl_max%kg
        !$OMP DO SCHEDULE(STATIC)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            cfl_tmp = abs((u(i,j,k)+u(i-1,j,k))*0.5_fp)*dxi + &
                      abs((v(i,j,k)+v(i,j-1,k))*0.5_fp)*dyi + &
                      abs((w(i,j,k)+w(i,j,k-1))*0.5_fp)*dzfi(k)
            if (cfl_tmp > cfl_max_thread) then
                cfl_max_thread = cfl_tmp
                i_tmp = i
                j_tmp = j
                k_tmp = k
            endif
        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP CRITICAL(find_cfl_max_among_threads)
        if (cfl_max_thread > cfl_max%value) then
            cfl_max%value = cfl_max_thread
            cfl_max%ig = i_tmp
            cfl_max%jg = j_tmp
            cfl_max%kg = k_tmp
        endif
        !$OMP END CRITICAL(find_cfl_max_among_threads)
        !$OMP END PARALLEL

        cfl_max%value = cfl_max%value*dt
        cfl_max%ig = cfl_max%ig + st(1) - 1
        cfl_max%jg = cfl_max%jg + st(2) - 1
        cfl_max%kg = cfl_max%kg + st(3) - 1

        cfl_max_root = cfl_max
        call MPI_REDUCE(cfl_max, cfl_max_root, 1, ValueIndexPair_type, MPI_MAX_INDEX, 0, &
                        MPI_COMM_WORLD, ierr)

        return
    end subroutine calcMaxCFL

    subroutine calcMaxDiv(nhalo, sz, dxi, dyi, dzfi, u, v, w, div_max_root)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        real(fp), intent(in) :: dxi, dyi
        real(fp), dimension(0:), intent(in) :: dzfi
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        type(value_index_pair_t), intent(out) :: div_max_root
        
        type(value_index_pair_t) :: div_max
        real(fp) :: div_tmp, div_max_thread
        integer :: i, j, k
        integer :: i_tmp, j_tmp, k_tmp

        div_max%value = 0.0_fp
        div_max%rank = myrank
        div_max%ig = 1
        div_max%jg = 1
        div_max%kg = 1
        !$OMP PARALLEL DEFAULT(SHARED) &
        !$OMP PRIVATE(i, j, k, div_tmp, div_max_thread, i_tmp, j_tmp, k_tmp)
        div_max_thread = div_max%value
        i_tmp = div_max%ig
        j_tmp = div_max%jg
        k_tmp = div_max%kg
        !$OMP DO SCHEDULE(STATIC)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            div_tmp = abs((u(i,j,k)-u(i-1,j,k))*dxi + &
                          (v(i,j,k)-v(i,j-1,k))*dyi + &
                          (w(i,j,k)-w(i,j,k-1))*dzfi(k))
            if (div_tmp > div_max_thread) then
                div_max_thread = div_tmp
                i_tmp = i
                j_tmp = j
                k_tmp = k
            endif
        enddo
        enddo
        enddo
        !$OMP END DO NOWAIT
        !$OMP CRITICAL(find_div_max_among_threads)
        if (div_max_thread > div_max%value) then
            div_max%value = div_max_thread
            div_max%ig = i_tmp
            div_max%jg = j_tmp
            div_max%kg = k_tmp
        endif
        !$OMP END CRITICAL(find_div_max_among_threads)
        !$OMP END PARALLEL

        div_max%ig = div_max%ig + st(1) - 1
        div_max%jg = div_max%jg + st(2) - 1
        div_max%kg = div_max%kg + st(3) - 1

        div_max_root = div_max
        call MPI_REDUCE(div_max, div_max_root, 1, ValueIndexPair_type, MPI_MAX_INDEX, 0, &
                        MPI_COMM_WORLD, ierr)

        return
    end subroutine calcMaxDiv

    subroutine checkCFL(cfl_limit, cfl_max, passed)
        implicit none
        real(fp), intent(in) :: cfl_limit
        type(value_index_pair_t), intent(in) :: cfl_max
        logical, intent(out) :: passed

        passed = .true.
        if (cfl_max%value > cfl_limit .or. cfl_max%value /= cfl_max%value) passed = .false.

        if ((myrank == 0) .and. (.not. passed)) then
            write(*,'(A)') 'PowerLLEL.ERROR.checkCFL: Maximum CFL number exceeds the threshold specified by the user!'
        endif

        return
    end subroutine checkCFL

    subroutine checkDiv(div_limit, div_max, passed)
        implicit none
        real(fp), intent(in) :: div_limit
        type(value_index_pair_t), intent(in) :: div_max
        logical, intent(out) :: passed

        passed = .true.
        if (div_max%value > div_limit .or. div_max%value /= div_max%value) passed = .false.
        
        if ((myrank == 0) .and. (.not. passed)) then
            write(*,'(A)') 'PowerLLEL.ERROR.checkDiv: Maximum divergence exceeds the threshold specified by the user!'
        endif

        return
    end subroutine checkDiv

    subroutine checkNaN(var, tag, passed_root)
        implicit none
        real(fp), dimension(:,:,:), intent(in) :: var
        character(*), intent(in), optional :: tag
        logical, intent(out) :: passed_root

        logical :: passed
        integer :: i, j, k
        
        passed = .true.

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = lbound(var,3), ubound(var,3)
        do j = lbound(var,2), ubound(var,2)
        !$  if (passed) then
                do i = lbound(var,1), ubound(var,1)
                    if (var(i,j,k) /= var(i,j,k)) then
        !$OMP           ATOMIC WRITE
                        passed = .false.
                        exit
                    endif
                enddo
        !$  else
        !$      exit
        !$  endif
        enddo
        enddo
        !$OMP END PARALLEL DO

        passed_root = passed
        call MPI_REDUCE(passed, passed_root, 1, MPI_LOGICAL, MPI_LAND, 0, MPI_COMM_WORLD, ierr)

        if ((myrank == 0) .and. (.not. passed_root)) then
            if (present(tag)) then
                write(*,'(A)') 'PowerLLEL.ERROR.checkNaN: NaN has been detected in <'//tag//'>!'
            else
                write(*,'(A)') 'PowerLLEL.ERROR.checkNaN: NaN has been detected!'
            endif
        endif

        return
    end subroutine checkNaN

end module mod_utils