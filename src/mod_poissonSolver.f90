module mod_poissonSolver
    use, intrinsic :: iso_c_binding
    use mod_type,  only: fp
    use mod_mpi,   only: neighbor_xyz, ierr, xsz, xst, ysz, yst, zsz, zst
#ifdef _PDD
    use mod_mpi,   only: MPI_REAL_FP, comm_cart_ypen, comm_cart_xpen
#endif
    use decomp_2d, only: transpose_x_to_y, transpose_y_to_x, &
                         transpose_y_to_z, transpose_z_to_y
    use mod_fft
    use mod_utils, only: abort
    !$ use omp_lib

#ifdef GPTL
    use gptl
#endif

    implicit none    

    include 'mpif.h'

    ! make everything private unless declared public
    private

    real(fp), dimension(0:1,3), save              :: rhs_cor
    real(fp), allocatable, dimension(:,:),   save :: lambdaxy
    real(fp), allocatable, dimension(:),     save :: a, c
#ifdef _PRECALC_TRID_COEFF
    real(fp), allocatable, dimension(:,:,:), save :: b
#else
    real(fp), allocatable, dimension(:),     save :: b
#endif
#ifdef _PDD
    real(fp), allocatable, dimension(:,:,:), save :: v_pdd, w_pdd
    real(fp), allocatable, dimension(:,:),   save :: y1_pdd, y2_pdd, y3_pdd
    real(fp), allocatable, dimension(:,:),   save :: tmp_var_pdd
    real(fp), allocatable, dimension(:,:),   save :: tmp_v_pdd
    integer,  save                                :: comm_cart_pdd
    integer,  save                                :: pdd_type
    integer                                       :: tag_pdd = 100
#endif
    real(fp), allocatable, dimension(:,:,:), target, save :: var_xpen, var_ypen, var_zpen
    integer, dimension(6), save                   :: neighbor_trid
#if defined(USE_MKL)
    type(fft_mkl_plan), dimension(2,2), save      :: fft_plan
#else
    type(c_ptr), dimension(2,2), save             :: fft_plan
#endif
    real(fp), save                                :: fft_normfactor
    logical,  save                                :: is_periodic_trid

#ifdef GPTL
    integer :: ret
#endif

    ! public user routines
    public :: initPoissonSolver, executePoissonSolver, freePoissonSolver

#ifdef USE_C
    interface
        subroutine init_poisson_solver(nx_global, ny_global, nz_global, dx, dy, dzf_global, &
                                       bctype_x, bctype_y, bctype_z, neighbor_xyz) &
            bind(C, name='init_poisson_solver')
            import
            integer(c_int), value :: nx_global, ny_global, nz_global
            real(c_double), value :: dx, dy
            real(c_double), dimension(*), intent(in) :: dzf_global
            character(kind=c_char), intent(in) :: bctype_x(*)
            character(kind=c_char), intent(in) :: bctype_y(*)
            character(kind=c_char), intent(in) :: bctype_z(*)
            integer(c_int), dimension(*), intent(in) :: neighbor_xyz
        end subroutine init_poisson_solver
        
        subroutine execute_poisson_solver(var) bind(C, name='execute_poisson_solver')
            import
            real(c_double), dimension(*) :: var
        end subroutine execute_poisson_solver

        subroutine free_poisson_solver() bind(C, name='free_poisson_solver')
            
        end subroutine free_poisson_solver
    end interface

    public :: init_poisson_solver, execute_poisson_solver, free_poisson_solver
#endif
    
contains
    
    !===================================================!
    !  Initialize basic data structures for the solver  !
    !===================================================!

    subroutine initPoissonSolver(nx_global, ny_global, nz_global, dx, dy, dzf_global, bctype, bcvalue)
        implicit none
        integer,  intent(in) :: nx_global
        integer,  intent(in) :: ny_global
        integer,  intent(in) :: nz_global
        real(fp), intent(in) :: dx
        real(fp), intent(in) :: dy
        real(fp), dimension(0:),    intent(in) :: dzf_global
        character(2), dimension(3), intent(in) :: bctype
        real(fp), dimension(0:1,3), intent(in) :: bcvalue
        
        integer,  dimension(3)   :: sz_trid, st_trid
        real(fp), dimension(0:1) :: dzc_b, dzf_b
        real(fp), allocatable, dimension(:) :: dzf
        real(fp), allocatable, dimension(:) :: lambdax, lambday
        real(fp), allocatable, dimension(:) :: b_tmp
        real(fp) :: a_tmp
        integer :: i, j, k, istatus

        ! allocate work arrays used in the transposition process
        allocate(var_xpen(xsz(1),xsz(2),xsz(3)), stat=istatus)
        if (istatus /= 0) call abort(105, "initPoissonSolver: Out of memory when allocating work arrays for transpositions!")
        allocate(var_ypen(ysz(1),ysz(2),ysz(3)), stat=istatus)
        if (istatus /= 0) call abort(105, "initPoissonSolver: Out of memory when allocating work arrays for transpositions!")
#ifndef _PDD
        allocate(var_zpen(zsz(1),zsz(2),zsz(3)), stat=istatus)
        if (istatus /= 0) call abort(105, "initPoissonSolver: Out of memory when allocating work arrays for transpositions!")
#endif

        ! determine a decomposition mode used for solving the tridiagonal system
#ifdef _PDD
        sz_trid = ysz; st_trid = yst
        neighbor_trid = neighbor_xyz(:,2)
#else
        sz_trid = zsz; st_trid = zst
        neighbor_trid = neighbor_xyz(:,3)
#endif
        
        ! initialize FFT
#if defined(USE_MKL)
        call initFFT(xsz, ysz, bctype(1:2), fft_plan, fft_normfactor)
#else
        call initFFT(xsz, ysz, bctype(1:2), var_xpen, var_ypen, fft_plan, fft_normfactor)
#endif
        
        ! calculate eigenvalues corresponding to BC types
        allocate(lambdax(sz_trid(1)))
        call getEigenvalues(st_trid(1), sz_trid(1), nx_global, bctype(1), lambdax)
        lambdax(:) = lambdax(:)/dx/dx
        allocate(lambday(sz_trid(2)))
        call getEigenvalues(st_trid(2), sz_trid(2), ny_global, bctype(2), lambday)
        lambday(:) = lambday(:)/dy/dy
        allocate(lambdaxy(sz_trid(1), sz_trid(2)))
        do j = 1, sz_trid(2)
        do i = 1, sz_trid(1)
            lambdaxy(i, j) = lambdax(i) + lambday(j)
        enddo
        enddo

        ! calculate coefficients of tridiagonal systems
#ifndef _PRECALC_TRID_COEFF
        allocate(b(sz_trid(3)), stat=istatus)
#else
        allocate(b(sz_trid(1), sz_trid(2), sz_trid(3)), stat=istatus)
        allocate(b_tmp(sz_trid(3)), stat=istatus)
#endif
        if (istatus /= 0) call abort(103, "initPoissonSolver: Out of memory when allocating coefficients of tridiagonal systems!")
        allocate(a(sz_trid(3)), stat=istatus)
        if (istatus /= 0) call abort(103, "initPoissonSolver: Out of memory when allocating coefficients of tridiagonal systems!")
        allocate(c(sz_trid(3)), stat=istatus)
        if (istatus /= 0) call abort(103, "initPoissonSolver: Out of memory when allocating coefficients of tridiagonal systems!")
        
        allocate(dzf(0:sz_trid(3)+1))
        dzf(:) = dzf_global(st_trid(3)-1:st_trid(3)+sz_trid(3))

#ifndef _PRECALC_TRID_COEFF
        call setTridCoeff(sz_trid(3), dzf, bctype(3), 'c', neighbor_trid(5:6), a, b, c)
#else
        call setTridCoeff(sz_trid(3), dzf, bctype(3), 'c', neighbor_trid(5:6), a, b_tmp, c)

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, a_tmp)
        do k = 1, sz_trid(3)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, sz_trid(2)
            do i = 1, sz_trid(1)
                b(i, j, k) = b_tmp(k) + lambdaxy(i, j)
            enddo
            enddo
            !$OMP END DO
        enddo

        ! decompose coefficient b
        do k = 2, sz_trid(3)
            !$OMP DO SCHEDULE(STATIC)
            do j = 1, sz_trid(2)
            do i = 1, sz_trid(1)
                a_tmp = a(k)/b(i, j, k-1)
                b(i, j, k) = b(i, j, k) - a_tmp*c(k-1)
            enddo
            enddo
            !$OMP END DO
        enddo
        !$OMP END PARALLEL
#endif

        ! determine whether the tridiagonal systems are periodic or not
        is_periodic_trid = .false.
        if (bctype(3) == 'PP') is_periodic_trid = .true.
        
        ! calculate the correction of the right-hand side according to BC types in x, y, z direction
        dzc_b(0) = 0.5_fp*(dzf_global(0)+dzf_global(1))
        dzc_b(1) = 0.5_fp*(dzf_global(nz_global)+dzf_global(nz_global+1))
        dzf_b(0) = dzf_global(1)
        dzf_b(1) = dzf_global(nz_global)
        call preprocessRHS(bctype, bcvalue, dx, dy, dzc_b, dzf_b)

#ifdef _PDD
        ! initialize MPI communication for PDD algorithm
        ! PDD comm in the bottom/top direction
        comm_cart_pdd = comm_cart_ypen
        call MPI_TYPE_VECTOR(sz_trid(1)*sz_trid(2), 1, 1, MPI_REAL_FP, pdd_type, ierr)
        call MPI_TYPE_COMMIT(pdd_type, ierr)
        
        ! initialize work arrays for PDD algorithm
        call initPDDArray(sz_trid)
#endif

        ! free temporary variables
        deallocate(lambdax)
        deallocate(lambday)
        deallocate(dzf)
#ifdef _PRECALC_TRID_COEFF
        deallocate(b_tmp)
#endif
        
        return
    end subroutine initPoissonSolver

#ifdef _PDD

    subroutine initPDDArray(sz)
        implicit none
        integer, dimension(3), intent(in) :: sz

        real(fp), dimension(sz(3)) :: b_tmp
        real(fp) :: a_tmp
        integer  :: i, j, k, ie, je, ke
        integer  :: istatus

        allocate(v_pdd (sz(1), sz(2), sz(3)), stat=istatus)
        if (istatus /= 0) call abort(104, "initPDDArray: Out of memory when allocating work arrays for PDD algorithm!")
        allocate(w_pdd (sz(1), sz(2), sz(3)), stat=istatus)
        if (istatus /= 0) call abort(104, "initPDDArray: Out of memory when allocating work arrays for PDD algorithm!")
        allocate(y1_pdd(sz(1), sz(2)), stat=istatus)
        if (istatus /= 0) call abort(104, "initPDDArray: Out of memory when allocating work arrays for PDD algorithm!")
        allocate(y2_pdd(sz(1), sz(2)), stat=istatus)
        if (istatus /= 0) call abort(104, "initPDDArray: Out of memory when allocating work arrays for PDD algorithm!")
        allocate(y3_pdd(sz(1), sz(2)), stat=istatus)
        if (istatus /= 0) call abort(104, "initPDDArray: Out of memory when allocating work arrays for PDD algorithm!")
        allocate(tmp_var_pdd(sz(1), sz(2)), stat=istatus)
        if (istatus /= 0) call abort(104, "initPDDArray: Out of memory when allocating work arrays for PDD algorithm!")
        allocate(tmp_v_pdd  (sz(1), sz(2)), stat=istatus)
        if (istatus /= 0) call abort(104, "initPDDArray: Out of memory when allocating work arrays for PDD algorithm!")

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,k)
        !$OMP DO SCHEDULE(STATIC)
        do k = 1, sz(3)
            v_pdd(:, :, k) = 0.0_fp
            w_pdd(:, :, k) = 0.0_fp
        enddo
        !$OMP END DO
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, sz(2)
            y1_pdd(:, j) = 0.0_fp
            y2_pdd(:, j) = 0.0_fp
            y3_pdd(:, j) = 0.0_fp
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        ie = sz(1); je = sz(2); ke = sz(3)

        do j = 1, je
        do i = 1, ie
            v_pdd(i, j, 1 ) = a(1 )
            w_pdd(i, j, ke) = c(ke)
        enddo
        enddo
    
        a(1 ) = 0.0_fp
        c(ke) = 0.0_fp

#ifdef _PRECALC_TRID_COEFF

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, a_tmp)
        do k = 2, ke
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            a_tmp = a(k)/b(i, j, k-1)
            v_pdd(i, j, k) = v_pdd(i, j, k) - a_tmp*v_pdd(i, j, k-1)
            w_pdd(i, j, k) = w_pdd(i, j, k) - a_tmp*w_pdd(i, j, k-1)
        enddo
        enddo
        !$OMP END DO
        enddo

        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            ! Important! To eliminate the singularity in the tridiagonal systems
            if (b(i, j, ke) /= 0.0) then
                v_pdd(i, j, ke) = v_pdd(i, j, ke)/b(i, j, ke)
                w_pdd(i, j, ke) = w_pdd(i, j, ke)/b(i, j, ke)
            else
                v_pdd(i, j, ke) = 0.0_fp
                w_pdd(i, j, ke) = 0.0_fp
            endif
        enddo
        enddo
        !$OMP END DO

        do k = ke-1, 1, -1
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            v_pdd(i, j, k) = (v_pdd(i, j, k) - c(k)*v_pdd(i, j, k+1))/b(i, j, k)
            w_pdd(i, j, k) = (w_pdd(i, j, k) - c(k)*w_pdd(i, j, k+1))/b(i, j, k)
        enddo
        enddo
        !$OMP END DO
        enddo
        !$OMP END PARALLEL

#else

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k, a_tmp, b_tmp)
        do j = 1, je
        do i = 1, ie
            b_tmp(:) = b(:) + lambdaxy(i, j)

            do k = 2, ke
                a_tmp = a(k)/b_tmp(k-1)
                b_tmp(k) = b_tmp(k) - a_tmp*c(k-1)
                v_pdd(i, j, k) = v_pdd(i, j, k) - a_tmp*v_pdd(i, j, k-1)
                w_pdd(i, j, k) = w_pdd(i, j, k) - a_tmp*w_pdd(i, j, k-1)
            enddo
            
            ! Important! To eliminate the singularity in the tridiagonal systems
            if (b_tmp(ke) /= 0.0) then
                v_pdd(i, j, ke) = v_pdd(i, j, ke)/b_tmp(ke)
                w_pdd(i, j, ke) = w_pdd(i, j, ke)/b_tmp(ke)
            else
                v_pdd(i, j, ke) = 0.0_fp
                w_pdd(i, j, ke) = 0.0_fp
            endif

            do k = ke-1, 1, -1
                v_pdd(i, j, k) = (v_pdd(i, j, k) - c(k)*v_pdd(i, j, k+1))/b_tmp(k)
                w_pdd(i, j, k) = (w_pdd(i, j, k) - c(k)*w_pdd(i, j, k+1))/b_tmp(k)
            enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

#endif
    
        call MPI_SENDRECV(  v_pdd(1,1,1), 1, pdd_type, neighbor_trid(5), tag_pdd, &
                          tmp_v_pdd(1,1), 1, pdd_type, neighbor_trid(6), tag_pdd, &
                          comm_cart_pdd, MPI_STATUS_IGNORE, ierr)
        
        return
    end subroutine initPDDArray

#endif

    subroutine setTridCoeff(n, dzf, bctype, c_or_f, neighbor, a, b, c)
        implicit none
        integer,  intent(in) :: n
        real(fp), dimension(0:), intent(in) :: dzf
        character(2),            intent(in) :: bctype
        character(1),            intent(in) :: c_or_f
        integer,  dimension(2),  intent(in) :: neighbor
        real(fp), dimension(n),  intent(out) :: a, b, c

        real(fp) :: factor
        integer :: k
        
        select case(c_or_f)
        case('c')
            do k = 1, n
                a(k) = 2.0_fp/( dzf(k)*(dzf(k-1)+dzf(k)) )
                c(k) = 2.0_fp/( dzf(k)*(dzf(k+1)+dzf(k)) )
            enddo
        case('f')
            do k = 1, n
                a(k) = 2.0_fp/( dzf(k)*(dzf(k+1)+dzf(k)) )
                c(k) = 2.0_fp/( dzf(k+1)*(dzf(k+1)+dzf(k)) )
            enddo
        end select

        b(:) = - a(:) - c(:)

        ! coefficients correction according to BC types
        if (neighbor(1) == MPI_PROC_NULL) then
            select case(bctype(1:1))
            case('P')
                factor = 0.0_fp
            case('D')
                factor = -1.0_fp
            case('N')
                factor = 1.0_fp
            end select

            select case(c_or_f)
            case('c')
                b(1) = b(1) + factor*a(1)
                a(1) = (abs(factor)-1.0_fp)*a(1)
            case('f')
                if(bctype(1:1) == 'N') then
                    b(1) = b(1) + factor*a(1)
                    a(1) = (abs(factor)-1.0_fp)*a(1)
                endif
            end select
        endif
        if (neighbor(2) == MPI_PROC_NULL) then
            select case(bctype(2:2))
            case('P')
                factor = 0.0_fp
            case('D')
                factor = -1.0_fp
            case('N')
                factor = 1.0_fp
            end select

            select case(c_or_f)
            case('c')
                b(n) = b(n) + factor*c(n)
                c(n) = (abs(factor)-1.0_fp)*c(n)
            case('f')
                if(bctype(2:2) == 'N') then
                    b(n) = b(n) + factor*c(n)
                    c(n) = (abs(factor)-1.0_fp)*c(n)
                endif
            end select
        endif

        return
    end subroutine setTridCoeff

    subroutine preprocessRHS(bctype, bcvalue, dx, dy, dzc_b, dzf_b)
        implicit none
        character(2), dimension(3), intent(in) :: bctype
        real(fp), dimension(0:1,3), intent(in) :: bcvalue
        real(fp),                   intent(in) :: dx, dy
        real(fp), dimension(0:1),   intent(in) :: dzc_b, dzf_b

        integer :: ibound, idir

        ! x direction
        idir = 1
        do ibound = 0, 1
            select case(bctype(idir)(ibound+1:ibound+1))
            case('P')
                rhs_cor(ibound, idir) = 0.0_fp
            case('D')
                rhs_cor(ibound, idir) = -2.0_fp*bcvalue(ibound, idir)/dx/dx
            case('N')
                if (ibound == 0) then
                    rhs_cor(ibound, idir) =  dx*bcvalue(ibound, idir)/dx/dx
                else
                    rhs_cor(ibound, idir) = -dx*bcvalue(ibound, idir)/dx/dx
                endif
            end select
        enddo
        ! y direction
        idir = 2
        do ibound = 0, 1
            select case(bctype(idir)(ibound+1:ibound+1))
            case('P')
                rhs_cor(ibound, idir) = 0.0_fp
            case('D')
                rhs_cor(ibound, idir) = -2.0_fp*bcvalue(ibound, idir)/dy/dy
            case('N')
                if (ibound == 0) then
                    rhs_cor(ibound, idir) =  dy*bcvalue(ibound, idir)/dy/dy
                else
                    rhs_cor(ibound, idir) = -dy*bcvalue(ibound, idir)/dy/dy
                endif
            end select
        enddo
        ! z direction
        idir = 3
        do ibound = 0, 1
            select case(bctype(idir)(ibound+1:ibound+1))
            case('P')
                rhs_cor(ibound, idir) = 0.0_fp
            case('D')
                rhs_cor(ibound, idir) = -2.0_fp*bcvalue(ibound, idir)/dzc_b(ibound)/dzf_b(ibound)
            case('N')
                if (ibound == 0) then
                    rhs_cor(ibound, idir) =  dzc_b(ibound)*bcvalue(ibound, idir)/dzc_b(ibound)/dzf_b(ibound)
                else
                    rhs_cor(ibound, idir) = -dzc_b(ibound)*bcvalue(ibound, idir)/dzc_b(ibound)/dzf_b(ibound)
                endif
            end select
        enddo

        return
    end subroutine preprocessRHS
    
    !=================================!
    !  Execute the solving procedure  !
    !=================================!
    
    subroutine executePoissonSolver(var)
        implicit none
        real(fp), dimension(0:xsz(1)+1,0:xsz(2)+1,0:xsz(3)+1), intent(inout) :: var
        
        real(fp), dimension(:,:,:), pointer :: y_ptr => null()
        real(fp), dimension(:,:,:), pointer :: z_ptr => null()
        integer :: i, j, k

#ifdef GPTL
        ret = gptlstart('--Copy input array')
#endif

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, xsz(3)
        do j = 1, xsz(2)
        do i = 1, xsz(1)
            var_xpen(i, j, k) = var(i, j, k)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

#ifdef GPTL
        ret = gptlstop('--Copy input array')
        ret = gptlstart('--Forward X-FFT')
#endif

        ! FFT along x direction (x-pencil); 
        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
        do k = 1, xsz(3)
        do j = 1, xsz(2)
            call executeFFT(fft_plan(1,1), var_xpen(:,j,k))
        enddo
        enddo

#ifdef GPTL
        ret = gptlstop('--Forward X-FFT')
        ret = gptlstart('--Transpose x to y')
#endif

        if (all(xsz == ysz)) then
            y_ptr => var_xpen
        else
            ! global transposition (x-pencil to y-pencil);
            call transpose_x_to_y(var_xpen, var_ypen)
            y_ptr => var_ypen
        endif

#ifdef GPTL
        ret = gptlstop('--Transpose x to y')
        ret = gptlstart('--Forward Y-FFT')
#endif

        ! FFT along y direction (y-pencil);
        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, k)
        do k = 1, ysz(3)
        do i = 1, ysz(1)
            call executeFFT(fft_plan(1,2), y_ptr(i,:,k))
        enddo
        enddo

#ifdef GPTL
        ret = gptlstop('--Forward Y-FFT')
#endif

#ifdef _PDD

#ifdef GPTL
        ret = gptlstart('--Solve trid')
#endif

        ! solve a series of tridiagonal systems in parallel (y-pencil);
        if (is_periodic_trid) then
            ! call pSolvePeriodicTrid(y_ptr)
        else
            call pSolveTrid(ysz, y_ptr)
        endif

#ifdef GPTL
        ret = gptlstop('--Solve trid')
#endif

#else

#ifdef GPTL
    ret = gptlstart('--Transpose y to z')
#endif

        if (all(ysz == zsz)) then
            z_ptr => y_ptr
        else
            ! global transposition (y-pencil to z-pencil);
            call transpose_y_to_z(y_ptr, var_zpen)
            z_ptr => var_zpen
        endif

#ifdef GPTL
        ret = gptlstop('--Transpose y to z')
        ret = gptlstart('--Solve trid')
#endif

        ! solve a series of tridiagonal systems in serial (z-pencil);
        if (is_periodic_trid) then
            ! call sSolvePeriodicTrid(z_ptr)
        else
            call sSolveTrid(zsz, z_ptr)
        endif

#ifdef GPTL
        ret = gptlstop('--Solve trid')
        ret = gptlstart('--Transpose z to y')
#endif

        if (any(ysz /= zsz)) then
            ! global transposition (z-pencil to y-pencil);
            call transpose_z_to_y(z_ptr, y_ptr)
        endif

#ifdef GPTL
        ret = gptlstop('--Transpose z to y')
#endif

#endif

#ifdef GPTL
        ret = gptlstart('--Backward Y-FFT')
#endif

        ! iFFT along y direction (y-pencil);
        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, k)
        do k = 1, ysz(3)
        do i = 1, ysz(1)
            call executeFFT(fft_plan(2,2), y_ptr(i,:,k))
        enddo
        enddo

#ifdef GPTL
        ret = gptlstop('--Backward Y-FFT')
        ret = gptlstart('--Transpose y to x')
#endif

        if (any(xsz /= ysz)) then
            ! global transposition (y-pencil to x-pencil);
            call transpose_y_to_x(y_ptr, var_xpen)
        endif

#ifdef GPTL
        ret = gptlstop('--Transpose y to x')
        ret = gptlstart('--Backward X-FFT')
#endif

        ! iFFT along x direction (x-pencil);
        !$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(j, k)
        do k = 1, xsz(3)
        do j = 1, xsz(2)
            call executeFFT(fft_plan(2,1), var_xpen(:,j,k))
        enddo
        enddo

#ifdef GPTL
        ret = gptlstop('--Backward X-FFT')
        ret = gptlstart('--Copy output array')
#endif

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, xsz(3)
        do j = 1, xsz(2)
        do i = 1, xsz(1)
            var(i, j, k) = var_xpen(i, j, k)*fft_normfactor
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

#ifdef GPTL
        ret = gptlstop('--Copy output array')
#endif

        nullify(y_ptr)
        nullify(z_ptr)
        
        return
    end subroutine executePoissonSolver

    subroutine correctRHS(sz, neighbor, rhs)
        implicit none
        integer,  dimension(3), intent(in) :: sz
        integer,  dimension(6), intent(in) :: neighbor
        real(fp), dimension(0:sz(1)+1,0:sz(2)+1,0:sz(3)+1), intent(inout) :: rhs

        integer :: i, j, k, ie, je, ke

        ie = sz(1); je = sz(2); ke = sz(3)

        ! x direction
        if (neighbor(1) == MPI_PROC_NULL) then
            do k = 1, ke
            do j = 1, je
                rhs(1 , j, k) = rhs(1 , j, k) + rhs_cor(0,1)
            enddo
            enddo
        endif
        if (neighbor(2) == MPI_PROC_NULL) then
            do k = 1, ke
            do j = 1, je
                rhs(ie, j, k) = rhs(ie, j, k) + rhs_cor(1,1)
            enddo
            enddo
        endif
        ! y direction
        if (neighbor(3) == MPI_PROC_NULL) then
            do k = 1, ke
            do i = 1, ie
                rhs(i, 1 , k) = rhs(i, 1 , k) + rhs_cor(0,2)
            enddo
            enddo
        endif
        if (neighbor(4) == MPI_PROC_NULL) then
            do k = 1, ke
            do i = 1, ie
                rhs(i, je, k) = rhs(i, je, k) + rhs_cor(1,2)
            enddo
            enddo
        endif
        ! z direction
        if (neighbor(5) == MPI_PROC_NULL) then
            do j = 1, je
            do i = 1, ie
                rhs(i, j, 1 ) = rhs(i, j, 1 ) + rhs_cor(0,3)
            enddo
            enddo
        endif
        if (neighbor(6) == MPI_PROC_NULL) then
            do j = 1, je
            do i = 1, ie
                rhs(i, j, ke) = rhs(i, j, ke) + rhs_cor(1,3)
            enddo
            enddo
        endif

        return
    end subroutine correctRHS

#ifdef _PDD

    subroutine pSolveTrid(sz, var)
        implicit none
        integer,  dimension(3), intent(in) :: sz
        real(fp), dimension(1:sz(1),1:sz(2),1:sz(3)), intent(inout) :: var
        
        real(fp), dimension(sz(3)) :: b_tmp
        real(fp) :: a_tmp, det_pdd
        integer  :: i, j, k, ie, je, ke

        ie = sz(1); je = sz(2); ke = sz(3)
        
#ifdef _PRECALC_TRID_COEFF

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, a_tmp)
        do k = 2, ke
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            a_tmp = a(k)/b(i, j, k-1)
            var(i, j, k) = var(i, j, k) - a_tmp*var(i, j, k-1)
        enddo
        enddo
        !$OMP END DO
        enddo

        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            ! Important! To eliminate the singularity in the tridiagonal systems
            if (b(i, j, ke) /= 0.0) then
                var(i, j, ke) = var(i, j, ke)/b(i, j, ke)
            else
                var(i, j, ke) = 0.0
            endif
        enddo
        enddo
        !$OMP END DO

        do k = ke-1, 1, -1
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            var(i, j, k) = (var(i, j, k) - c(k)*var(i, j, k+1))/b(i, j, k)
        enddo
        enddo
        !$OMP END DO
        enddo
        !$OMP END PARALLEL

#else

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k, a_tmp, b_tmp)
        do j = 1, je
        do i = 1, ie
            b_tmp(:) = b(:) + lambdaxy(i, j)

            do k = 2, ke
                a_tmp = a(k)/b_tmp(k-1)
                b_tmp(k) = b_tmp(k) - a_tmp*c(k-1)
                var(i, j, k) = var(i, j, k) - a_tmp*var(i, j, k-1)
            enddo
            
            ! Important! To eliminate the singularity in the tridiagonal systems
            if (b_tmp(ke) /= 0.0) then
                var(i, j, ke) = var(i, j, ke)/b_tmp(ke)
            else
                var(i, j, ke) = 0.0_fp
            endif

            do k = ke-1, 1, -1
                var(i, j, k) = (var(i, j, k) - c(k)*var(i, j, k+1))/b_tmp(k)
            enddo

        enddo
        enddo
        !$OMP END PARALLEL DO

#endif

#ifdef GPTL
        ret = gptlstart('----Comm in PDD')
#endif

        call MPI_SENDRECV(      var(1,1,1), 1, pdd_type, neighbor_trid(5), tag_pdd, &
                          tmp_var_pdd(1,1), 1, pdd_type, neighbor_trid(6), tag_pdd, &
                          comm_cart_pdd, MPI_STATUS_IGNORE, ierr)
        
#ifdef GPTL
        ret = gptlstop('----Comm in PDD')
#endif

        if (neighbor_trid(6) /= MPI_PROC_NULL) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k, det_pdd)
            do j = 1, je
            do i = 1, ie
                det_pdd = w_pdd(i, j, ke)*tmp_v_pdd(i, j) - 1.0_fp
                y2_pdd(i, j) = (var(i, j, ke)*tmp_v_pdd(i, j) - tmp_var_pdd(i, j))/det_pdd
                y3_pdd(i, j) = (tmp_var_pdd(i, j)*w_pdd(i, j, ke) - var(i, j, ke))/det_pdd
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif

#ifdef GPTL
        ret = gptlstart('----Comm in PDD')
#endif

        call MPI_SENDRECV(y3_pdd(1,1), 1, pdd_type, neighbor_trid(6), tag_pdd, &
                          y1_pdd(1,1), 1, pdd_type, neighbor_trid(5), tag_pdd, &
                          comm_cart_pdd, MPI_STATUS_IGNORE, ierr)

#ifdef GPTL
        ret = gptlstop('----Comm in PDD')
#endif

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, ke
        do j = 1, je
        do i = 1, ie
            var(i, j, k) = var(i, j, k) - v_pdd(i, j, k)*y1_pdd(i, j) - w_pdd(i, j, k)*y2_pdd(i, j)
        enddo
        enddo
        enddo

        return
    end subroutine pSolveTrid

#endif

    subroutine sSolveTrid(sz, var)
        implicit none
        integer,  dimension(3), intent(in) :: sz
        real(fp), dimension(1:sz(1),1:sz(2),1:sz(3)), intent(inout) :: var

        real(fp), dimension(sz(3)) :: b_tmp
        real(fp) :: a_tmp
        integer  :: i, j, k, ie, je, ke

        ie = sz(1); je = sz(2); ke = sz(3)
        
#ifdef _PRECALC_TRID_COEFF

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, a_tmp)
        do k = 2, ke
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            a_tmp = a(k)/b(i, j, k-1)
            var(i, j, k) = var(i, j, k) - a_tmp*var(i, j, k-1)
        enddo
        enddo
        !$OMP END DO
        enddo

        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            ! Important! To eliminate the singularity in the tridiagonal systems
            if (b(i, j, ke) /= 0.0) then
                var(i, j, ke) = var(i, j, ke)/b(i, j, ke)
            else
                var(i, j, ke) = 0.0
            endif
        enddo
        enddo
        !$OMP END DO

        do k = ke-1, 1, -1
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, je
        do i = 1, ie
            var(i, j, k) = (var(i, j, k) - c(k)*var(i, j, k+1))/b(i, j, k)
        enddo
        enddo
        !$OMP END DO
        enddo
        !$OMP END PARALLEL

#else

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k, a_tmp, b_tmp)
        do j = 1, je
        do i = 1, ie
            b_tmp(:) = b(:) + lambdaxy(i, j)

            do k = 2, ke
                a_tmp = a(k)/b_tmp(k-1)
                b_tmp(k) = b_tmp(k) - a_tmp*c(k-1)
                var(i, j, k) = var(i, j, k) - a_tmp*var(i, j, k-1)
            enddo
            
            ! Important! To eliminate the singularity in the tridiagonal systems
            if (b_tmp(ke) /= 0.0) then
                var(i, j, ke) = var(i, j, ke)/b_tmp(ke)
            else
                var(i, j, ke) = 0.0_fp
            endif

            do k = ke-1, 1, -1
                var(i, j, k) = (var(i, j, k) - c(k)*var(i, j, k+1))/b_tmp(k)
            enddo

        enddo
        enddo
        !$OMP END PARALLEL DO

#endif

        return
    end subroutine sSolveTrid
    
    !=========================================!
    !  Release the memory used by the solver  !
    !=========================================!
    
    subroutine freePoissonSolver()
        implicit none

        ! release work arrays
        if (allocated(var_xpen)) deallocate(var_xpen)
        if (allocated(var_ypen)) deallocate(var_ypen)
        if (allocated(var_zpen)) deallocate(var_zpen)
        
        ! release tridiagonal coefficients arrays
        deallocate(lambdaxy)
        deallocate(a, b, c)

#ifdef _PDD
        ! release MPI objects
        call MPI_TYPE_FREE(pdd_type, ierr)

        ! release PDD related arrays
        deallocate(v_pdd, w_pdd)
        deallocate(y1_pdd, y2_pdd, y3_pdd)
        deallocate(tmp_var_pdd)
        deallocate(tmp_v_pdd)
#endif

        ! release FFT
        call freeFFT(fft_plan)

        return
    end subroutine freePoissonSolver

    
end module mod_poissonSolver