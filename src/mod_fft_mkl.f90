include "mkl_dfti.f90"
include "mkl_trig_transforms.f90"

#define UNROLL_F 4
!#define FFT_MKL_CESTOR_CC

module fft_mkl
    use, intrinsic :: iso_c_binding
    use mkl_dfti
    use mkl_trig_transforms
    implicit none
    ! make everything private unless declared public
    private

    ! FFTW KIND constants borrowed from the header file fftw3.f03
    integer(c_int), parameter, public :: FFTW_R2HC = 0
    integer(c_int), parameter, public :: FFTW_HC2R = 1
    integer(c_int), parameter, public :: FFTW_DHT = 2
    integer(c_int), parameter, public :: FFTW_REDFT00 = 3
    integer(c_int), parameter, public :: FFTW_REDFT01 = 4
    integer(c_int), parameter, public :: FFTW_REDFT10 = 5
    integer(c_int), parameter, public :: FFTW_REDFT11 = 6
    integer(c_int), parameter, public :: FFTW_RODFT00 = 7
    integer(c_int), parameter, public :: FFTW_RODFT01 = 8
    integer(c_int), parameter, public :: FFTW_RODFT10 = 9
    integer(c_int), parameter, public :: FFTW_RODFT11 = 10
    integer(c_int), parameter, public :: FFTW_FORWARD = -1
    integer(c_int), parameter, public :: FFTW_BACKWARD = +1

    integer, parameter, public :: MKL_PACKED_FMT_PACK  = 11
    integer, parameter, public :: MKL_PACKED_FMT_PERM  = 12
    integer, parameter, public :: MKL_PACKED_FMT_PERM2 = 13

    type fft_mkl_plan_desc
        type(DFTI_DESCRIPTOR), pointer :: dfti_desc => null()
        integer :: fft_type   ! 1 - dfti in DFTI_COMPLEX_REAL storage; 2 - dfti in DFTI_COMPLEX_COMPLEX storage; 3 - tt
        integer :: packed_fmt ! 11 - DFTI_PACK_FORMAT; 12 - DFTI_PERM_FORMAT; 13 - PERM2 format (new); -1 - unused
        integer :: n, howmany
        integer :: istride, ostride
        integer :: idist, odist
        integer,  allocatable, dimension(:) :: ipar
        real(c_double), allocatable, dimension(:) :: dpar
        real(c_double), allocatable, dimension(:) :: buffer
    end type

    type fft_mkl_plan
        type(fft_mkl_plan_desc) :: desc
        procedure(prepost_interface), pointer, nopass :: preproc => null()
        procedure(execute_interface), pointer, nopass :: execute => null()
        procedure(prepost_interface), pointer, nopass :: posproc => null()
        procedure(destroy_interface), pointer, nopass :: destroy => null()
    end type

    abstract interface
        subroutine prepost_interface(desc, work)
            import :: fft_mkl_plan_desc, c_double
            type(fft_mkl_plan_desc), intent(inout) :: desc
            real(c_double), dimension(*), intent(inout) :: work
        end subroutine prepost_interface

        subroutine execute_interface(desc, work)
            import :: fft_mkl_plan_desc, c_double
            type(fft_mkl_plan_desc), intent(inout) :: desc
            real(c_double), dimension(*), intent(inout) :: work
        end subroutine execute_interface

        subroutine destroy_interface(desc)
            import :: fft_mkl_plan_desc
            type(fft_mkl_plan_desc), intent(inout) :: desc
        end subroutine destroy_interface
    end interface

    public :: fft_mkl_plan, fft_mkl_plan_many_r2r_1d, fft_mkl_execute_r2r, fft_mkl_destroy_plan

contains

    subroutine fft_mkl_plan_many_r2r_1d(plan, packed_fmt, fftw_kind, n, howmany, istride, idist, ostride, odist)
        implicit none
        type(fft_mkl_plan), intent(out) :: plan
        integer, intent(in) :: packed_fmt
        integer(c_int), intent(in) :: fftw_kind
        integer, intent(in) :: n, howmany
        integer, intent(in) :: istride, idist
        integer, intent(in) :: ostride, odist

        ! Workarounds, not equivalent to FFTW, but useful for solving PDE by eigenfunction expansion method.
        ! Note that the vector containing corresponding eigen values should be reorganized to match with the MKL packed format.

        select case(fftw_kind)
        case(FFTW_DHT)
#ifdef FFT_MKL_CESTOR_CC
            call dfti_plan_cestorage_cc(plan, FFTW_FORWARD, packed_fmt, n, howmany, istride, idist, ostride, odist)
#else
            call dfti_plan_cestorage_cr(plan, FFTW_FORWARD, packed_fmt, n, howmany, istride, idist, ostride, odist)
#endif
        case(FFTW_R2HC)
#ifdef FFT_MKL_CESTOR_CC
            call dfti_plan_cestorage_cc(plan, FFTW_FORWARD, packed_fmt, n, howmany, istride, idist, ostride, odist)
#else
            call dfti_plan_cestorage_cr(plan, FFTW_FORWARD, packed_fmt, n, howmany, istride, idist, ostride, odist)
#endif
        case(FFTW_HC2R)
#ifdef FFT_MKL_CESTOR_CC
            call dfti_plan_cestorage_cc(plan, FFTW_BACKWARD, packed_fmt, n, howmany, istride, idist, ostride, odist)
#else
            call dfti_plan_cestorage_cr(plan, FFTW_BACKWARD, packed_fmt, n, howmany, istride, idist, ostride, odist)
#endif
        case(FFTW_REDFT10)
            call tt_plan(plan, MKL_STAGGERED_COSINE_TRANSFORM, FFTW_BACKWARD, n, howmany, istride, idist, ostride, odist)
        case(FFTW_REDFT01)
            call tt_plan(plan, MKL_STAGGERED_COSINE_TRANSFORM, FFTW_FORWARD, n, howmany, istride, idist, ostride, odist)
        case(FFTW_REDFT11)
            call tt_plan(plan, MKL_STAGGERED2_COSINE_TRANSFORM, FFTW_BACKWARD, n, howmany, istride, idist, ostride, odist)
        case(FFTW_RODFT10)
            call tt_plan(plan, MKL_STAGGERED_SINE_TRANSFORM, FFTW_BACKWARD, n, howmany, istride, idist, ostride, odist)
        case(FFTW_RODFT01)
            call tt_plan(plan, MKL_STAGGERED_SINE_TRANSFORM, FFTW_FORWARD, n, howmany, istride, idist, ostride, odist)
        case(FFTW_RODFT11)
            call tt_plan(plan, MKL_STAGGERED2_SINE_TRANSFORM, FFTW_BACKWARD, n, howmany, istride, idist, ostride, odist)
        end select

        return
    end subroutine fft_mkl_plan_many_r2r_1d
    
    subroutine dfti_plan_cestorage_cc(plan, fwd_or_bwd, packed_fmt, n, howmany, istride, idist, ostride, odist)
        implicit none
        type(fft_mkl_plan), intent(out) :: plan
        integer, intent(in) :: fwd_or_bwd
        integer, intent(in) :: packed_fmt
        integer, intent(in) :: n, howmany
        integer, intent(in) :: istride, idist
        integer, intent(in) :: ostride, odist

        integer :: strides(2), status

        plan%desc%packed_fmt = packed_fmt
        plan%desc%fft_type = 2
        plan%desc%n       = n
        plan%desc%howmany = howmany
        plan%desc%istride = istride
        plan%desc%ostride = ostride
        plan%desc%idist = idist
        plan%desc%odist = odist
        allocate(plan%desc%buffer( (n+2)*UNROLL_F ))
        
        status = DftiCreateDescriptor(plan%desc%dfti_desc, DFTI_DOUBLE, DFTI_REAL, 1, n)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_NUMBER_OF_TRANSFORMS, 1)
        strides(1) = 0
        strides(2) = 1
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_INPUT_STRIDES,  strides)
        strides(2) = 1
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_OUTPUT_STRIDES, strides)
        status = DftiCommitDescriptor(plan%desc%dfti_desc)
        
        if (fwd_or_bwd == FFTW_FORWARD) then
            
            plan%execute => dfti_execute_fi

            select case(packed_fmt)
            case(MKL_PACKED_FMT_PACK)
                plan%posproc => dfti_fwd_posproc_pack
            case(MKL_PACKED_FMT_PERM)
                plan%posproc => dfti_fwd_posproc_perm
            case(MKL_PACKED_FMT_PERM2)
                plan%posproc => dfti_fwd_posproc_perm2
            end select

        else if (fwd_or_bwd == FFTW_BACKWARD) then
            
            plan%execute => dfti_execute_bi
            
            select case(packed_fmt)
            case(MKL_PACKED_FMT_PACK)
                plan%preproc => dfti_bwd_preproc_pack
            case(MKL_PACKED_FMT_PERM)
                plan%preproc => dfti_bwd_preproc_perm
            case(MKL_PACKED_FMT_PERM2)
                plan%preproc => dfti_bwd_preproc_perm2
            end select

        endif

        plan%destroy => dfti_destroy

        return
    end subroutine dfti_plan_cestorage_cc
    
    subroutine dfti_plan_cestorage_cr(plan, fwd_or_bwd, packed_fmt, n, howmany, istride, idist, ostride, odist)
        implicit none
        type(fft_mkl_plan), intent(out) :: plan
        integer, intent(in) :: fwd_or_bwd
        integer, intent(in) :: packed_fmt
        integer, intent(in) :: n, howmany
        integer, intent(in) :: istride, idist
        integer, intent(in) :: ostride, odist

        integer :: strides(2), status, fmt

        select case(packed_fmt)
        case(MKL_PACKED_FMT_PACK)
            plan%desc%packed_fmt = packed_fmt
            fmt = DFTI_PACK_FORMAT
        case(MKL_PACKED_FMT_PERM)
            plan%desc%packed_fmt = packed_fmt
            fmt = DFTI_PERM_FORMAT
        case(MKL_PACKED_FMT_PERM2)
            plan%desc%packed_fmt = MKL_PACKED_FMT_PERM ! MKL does not support PERM2 format
            fmt = DFTI_PERM_FORMAT
        end select
        plan%desc%fft_type = 1
        plan%desc%n       = n
        plan%desc%howmany = howmany
        plan%desc%istride = istride
        plan%desc%ostride = ostride
        plan%desc%idist = idist
        plan%desc%odist = odist
        
        status = DftiCreateDescriptor(plan%desc%dfti_desc, DFTI_DOUBLE, DFTI_REAL, 1, n)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_PACKED_FORMAT, fmt)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_NUMBER_OF_TRANSFORMS, howmany)
        strides(1) = 0
        strides(2) = istride
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_INPUT_STRIDES,  strides)
        strides(2) = ostride
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_OUTPUT_STRIDES, strides)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_INPUT_DISTANCE,  idist)
        status = DftiSetValue(plan%desc%dfti_desc, DFTI_OUTPUT_DISTANCE, odist)
        status = DftiCommitDescriptor(plan%desc%dfti_desc)
        
        if (fwd_or_bwd == FFTW_FORWARD) then
            plan%execute => dfti_execute_fi
        else if (fwd_or_bwd == FFTW_BACKWARD) then
            plan%execute => dfti_execute_bi
        endif

        plan%destroy => dfti_destroy

        return
    end subroutine dfti_plan_cestorage_cr

    subroutine dfti_fwd_posproc_pack(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: i, ist, ien

        ist = 2
        ien = desc%n
        work(ist : ien) = work(ist+1 : ien+1)

        return
    end subroutine dfti_fwd_posproc_pack

    subroutine dfti_fwd_posproc_perm(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: i, ist, ien

        ist = 2
        ien = desc%n
        if (mod(desc%n, 2) == 0) then
            work(ist) = work(ien+1)
        else
            work(ist : ien) = work(ist+1 : ien+1)
        endif

        return
    end subroutine dfti_fwd_posproc_perm

    subroutine dfti_fwd_posproc_perm2(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: i, ist, ien

        ist = 2
        ien = desc%n
        work(ist) = work(ien+1)

        return
    end subroutine dfti_fwd_posproc_perm2

    subroutine dfti_bwd_preproc_pack(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: ist, ien

        ist = 2
        ien = desc%n

        work(ist+1 : ien+1) = work(ist : ien)
        work(ist) = 0.0
        if (mod(desc%n,2) == 0) work(desc%n+2) = 0.0
        
        return
    end subroutine dfti_bwd_preproc_pack

    subroutine dfti_bwd_preproc_perm(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: ist, ien

        ist = 2
        ien = desc%n
        if (mod(desc%n,2) == 0) then
            work(ien+1) = work(ist)
            work(ien+2) = 0.0
        else
            work(ist+1 : ien+1) = work(ist : ien)
        endif
        work(ist) = 0.0
        
        return
    end subroutine dfti_bwd_preproc_perm

    subroutine dfti_bwd_preproc_perm2(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: ist, ien

        ist = 2
        ien = desc%n
        work(ien+1) = work(ist)
        work(ien+2) = 0.0
        work(ist) = 0.0
        
        return
    end subroutine dfti_bwd_preproc_perm2

    subroutine dfti_execute_fi(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: status

        status = DftiComputeForward(desc%dfti_desc, work)

        return
    end subroutine dfti_execute_fi

    subroutine dfti_execute_bi(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: status

        status = DftiComputeBackward(desc%dfti_desc, work)

        return
    end subroutine dfti_execute_bi

    subroutine dfti_destroy(desc)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        integer :: status

        status = DftiFreeDescriptor(desc%dfti_desc)

        return
    end subroutine dfti_destroy

    subroutine tt_plan(plan, tt_type, fwd_or_bwd, n, howmany, istride, idist, ostride, odist)
        implicit none
        type(fft_mkl_plan), intent(out) :: plan
        integer, intent(in) :: tt_type
        integer, intent(in) :: fwd_or_bwd
        integer, intent(in) :: n, howmany
        integer, intent(in) :: istride, idist
        integer, intent(in) :: ostride, odist

        real(c_double), allocatable, dimension(:) :: f
        integer :: ir

        plan%desc%fft_type = 3
        plan%desc%packed_fmt = -1 ! unused
        plan%desc%n       = n
        plan%desc%howmany = howmany
        plan%desc%istride = istride
        plan%desc%ostride = ostride
        plan%desc%idist = idist
        plan%desc%odist = odist

        allocate(plan%desc%ipar(128))
        allocate(plan%desc%dpar(5*n/2+2))
        call d_init_trig_transform(n, tt_type, plan%desc%ipar, plan%desc%dpar, ir)
        plan%desc%ipar(11) = 1  ! compatibility with FFTW, unnormalized
        allocate(f(n))
        call d_commit_trig_transform(f, plan%desc%dfti_desc, plan%desc%ipar, plan%desc%dpar, ir)
        deallocate(f)

        if ( istride > 1 ) then
            allocate(plan%desc%buffer(UNROLL_F*n))
        end if
        ! plan%preproc => tt_preproc
        ! plan%posproc => tt_posproc

        if (fwd_or_bwd == FFTW_FORWARD) then
            plan%execute => tt_execute_f
        else if (fwd_or_bwd == FFTW_BACKWARD) then
            plan%execute => tt_execute_b
        endif

        plan%destroy => tt_destroy

        return
    end subroutine tt_plan

    subroutine tt_preproc(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work

        return
    end subroutine tt_preproc

    subroutine tt_posproc(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work

        return
    end subroutine tt_posproc

    subroutine tt_execute_f(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: ir

        call d_forward_trig_transform(work, desc%dfti_desc, desc%ipar, desc%dpar, ir)

        return
    end subroutine tt_execute_f

    subroutine tt_execute_b(desc, work)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        real(c_double), dimension(*), intent(inout) :: work
        integer :: ir

        call d_backward_trig_transform(work, desc%dfti_desc, desc%ipar, desc%dpar, ir)

        return
    end subroutine tt_execute_b
    
    subroutine tt_destroy(desc)
        implicit none
        type(fft_mkl_plan_desc), intent(inout) :: desc
        integer :: ir

        call free_trig_transform(desc%dfti_desc, desc%ipar, ir)
        if (allocated(desc%ipar)) deallocate(desc%ipar)
        if (allocated(desc%dpar)) deallocate(desc%dpar)

        return
    end subroutine tt_destroy

    subroutine fft_mkl_execute_r2r(plan, work)
        implicit none
        type(fft_mkl_plan), intent(inout) :: plan
        real(c_double), dimension(*), intent(inout) :: work

        integer :: pos, i, j, irest, buf_lead_dim, iloc

        ! execute in-place FFT
        select case(plan%desc%fft_type)
        case(1)
            ! dfti in DFTI_COMPLEX_REAL storage

            call plan%execute(plan%desc, work)

        case(2)
            ! dfti in DFTI_COMPLEX_COMPLEX storage

            buf_lead_dim = plan%desc%n + 2
            irest = mod(plan%desc%howmany, UNROLL_F)
            i = 1
            do while(i <= irest)
                pos = (i-1) * plan%desc%idist + 1
                i = i + 1
                call gather(1, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos), buf_lead_dim, plan%desc%buffer)
                if(associated(plan%preproc)) call plan%preproc(plan%desc, plan%desc%buffer)
                call plan%execute(plan%desc, plan%desc%buffer)
                if(associated(plan%posproc)) call plan%posproc(plan%desc, plan%desc%buffer)
                call scatter(1, buf_lead_dim, plan%desc%buffer, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos))
            enddo
            do while(i <= plan%desc%howmany)
                pos = (i-1) * plan%desc%idist + 1
                i = i + UNROLL_F
                call gather(UNROLL_F, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos), buf_lead_dim, plan%desc%buffer)
                !GCC$ UNROLL UNROLL_F
                !DIR$ UNROLL(UNROLL_F)
                do j = 0, UNROLL_F-1
                    iloc = 1 + buf_lead_dim * j
                    if(associated(plan%preproc)) call plan%preproc(plan%desc, plan%desc%buffer(iloc))
                    call plan%execute(plan%desc, plan%desc%buffer(iloc))
                    if(associated(plan%posproc)) call plan%posproc(plan%desc, plan%desc%buffer(iloc))
                enddo
                call scatter(UNROLL_F, buf_lead_dim, plan%desc%buffer, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos))
            enddo
            ! No-unroll version
            ! buf_lead_dim = plan%desc%n + 2
            ! do i = 1, plan%desc%howmany
            !     pos = (i-1) * plan%desc%idist + 1
            !     call gather(1, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos), buf_lead_dim, plan%desc%buffer)
            !     if(associated(plan%preproc)) call plan%preproc(plan%desc, plan%desc%buffer)
            !     call plan%execute(plan%desc, plan%desc%buffer)
            !     if(associated(plan%posproc)) call plan%posproc(plan%desc, plan%desc%buffer)
            !     call scatter(1, buf_lead_dim, plan%desc%buffer, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos))
            ! enddo
            
        case(3)
            ! tt
            
            if ( plan%desc%istride > 1 ) then
                buf_lead_dim = plan%desc%n
                irest = mod(plan%desc%howmany, UNROLL_F)
                i = 1
                do while(i <= irest)
                    pos = (i-1) * plan%desc%idist + 1
                    i = i + 1
                    call gather(1, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos), buf_lead_dim, plan%desc%buffer)
                    if(associated(plan%preproc)) call plan%preproc(plan%desc, plan%desc%buffer)
                    call plan%execute(plan%desc, plan%desc%buffer)
                    if(associated(plan%posproc)) call plan%posproc(plan%desc, plan%desc%buffer)
                    call scatter(1, buf_lead_dim, plan%desc%buffer, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos))
                enddo
                do while(i <= plan%desc%howmany)
                    pos = (i-1) * plan%desc%idist + 1
                    i = i + UNROLL_F
                    call gather(UNROLL_F, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos), buf_lead_dim, plan%desc%buffer)
                    !GCC$ UNROLL UNROLL_F
                    !DIR$ UNROLL(UNROLL_F)
                    do j = 0, UNROLL_F-1
                        iloc = 1 + buf_lead_dim * j
                        if(associated(plan%preproc)) call plan%preproc(plan%desc, plan%desc%buffer(iloc))
                        call plan%execute(plan%desc, plan%desc%buffer(iloc))
                        if(associated(plan%posproc)) call plan%posproc(plan%desc, plan%desc%buffer(iloc))
                    enddo
                    call scatter(UNROLL_F, buf_lead_dim, plan%desc%buffer, plan%desc%n, plan%desc%istride, plan%desc%idist, work(pos))
                enddo
            else
                do i = 1, plan%desc%howmany
                    pos = (i-1) * plan%desc%idist + 1
                    if(associated(plan%preproc)) call plan%preproc(plan%desc, work(pos))
                    call plan%execute(plan%desc, work(pos))
                    if(associated(plan%posproc)) call plan%posproc(plan%desc, work(pos))
                enddo
            end if
        end select

        return
    end subroutine fft_mkl_execute_r2r

    subroutine gather(unroll_factor, n, stride, dist, src, lead_dim, dst)
        implicit none
        integer, intent(in) :: unroll_factor, n, stride, dist
        real(c_double), dimension(*), intent(in)  :: src
        integer, intent(in) :: lead_dim
        real(c_double), dimension(*), intent(out) :: dst

        integer :: i, j, pos

        if (stride == 1) then
            do j = 0, unroll_factor-1
                do i = 1, n
                   dst(i + lead_dim*j) = src(i + dist*j)
                end do
            end do
        else
            pos = 1
            do i = 1, n
                !GCC$ UNROLL UNROLL_F
                !DIR$ UNROLL(UNROLL_F)
                do j = 0, unroll_factor-1
                    dst(i + lead_dim*j) = src(pos + dist*j)
                end do
                pos = pos + stride
            end do
        endif

        return
    end subroutine gather

    subroutine scatter(unroll_factor, lead_dim, src, n, stride, dist, dst)
        implicit none
        integer, intent(in) :: lead_dim
        real(c_double), dimension(*), intent(in)  :: src
        integer, intent(in) :: unroll_factor, n, stride, dist
        real(c_double), dimension(*), intent(out) :: dst

        integer :: i, j, pos

        if (stride == 1) then
            do j = 0, unroll_factor-1
                do i = 1, n
                    dst(i + dist*j) = src(i + lead_dim*j)
                end do
            end do
        else
            pos = 1
            do i = 1, n
                !GCC$ UNROLL UNROLL_F
                !DIR$ UNROLL(UNROLL_F)
                do j = 0, unroll_factor-1
                    dst(pos + dist*j) = src(i + lead_dim*j)
                end do
                pos = pos + stride
            end do
        endif

        return
    end subroutine scatter

    subroutine fft_mkl_destroy_plan(plan)
        implicit none
        type(fft_mkl_plan), intent(inout) :: plan

        call plan%destroy(plan%desc)

        return
    end subroutine fft_mkl_destroy_plan
    
end module fft_mkl

module mod_fft
    use, intrinsic :: iso_c_binding
    use mod_type, only: fp
    !$ use omp_lib
    use fft_mkl
    
    implicit none
    private

    integer, parameter :: mkl_packed_fmt = MKL_PACKED_FMT_PERM

    public :: fft_mkl_plan, initFFT, getEigenvalues, executeFFT, freeFFT

contains
    subroutine initFFT(xsz, ysz, bctype_xy, fft_plan, fft_normfactor)
        implicit none
        integer, dimension(3), intent(in) :: xsz, ysz
        character(2), dimension(2), intent(in) :: bctype_xy
        type(fft_mkl_plan), dimension(2,2), intent(out) :: fft_plan
        real(fp), intent(out) :: fft_normfactor

        integer(c_int) :: kind_fwd, kind_bwd
        real(fp) :: normfactor_x, normfactor_y

        ! init single in-place FFT along x direction (x-pencil), idist & odist is unused here
        call getFFTWKind(bctype_xy(1), kind_fwd, kind_bwd)
        call getNormFactor(bctype_xy(1), normfactor_x)
        call fft_mkl_plan_many_r2r_1d(fft_plan(1,1), mkl_packed_fmt, kind_fwd, xsz(1), 1, 1, xsz(1), 1, xsz(1))
        call fft_mkl_plan_many_r2r_1d(fft_plan(2,1), mkl_packed_fmt, kind_bwd, xsz(1), 1, 1, xsz(1), 1, xsz(1))
        
        ! init single in-place FFT along y direction (y-pencil), idist & odist is unused here
        call getFFTWKind(bctype_xy(2), kind_fwd, kind_bwd)
        call getNormFactor(bctype_xy(2), normfactor_y)
        ! call fft_mkl_plan_many_r2r_1d(fft_plan(1,2), mkl_packed_fmt, kind_fwd, ysz(2), 1, ysz(1), 1, ysz(1), 1)
        ! call fft_mkl_plan_many_r2r_1d(fft_plan(2,2), mkl_packed_fmt, kind_bwd, ysz(2), 1, ysz(1), 1, ysz(1), 1)
        call fft_mkl_plan_many_r2r_1d(fft_plan(1,2), mkl_packed_fmt, kind_fwd, ysz(2), 1, 1, ysz(1), 1, ysz(1))
        call fft_mkl_plan_many_r2r_1d(fft_plan(2,2), mkl_packed_fmt, kind_bwd, ysz(2), 1, 1, ysz(1), 1, ysz(1))

        fft_normfactor = 1.0_fp/(normfactor_x*xsz(1)*normfactor_y*ysz(2))

        return
    end subroutine initFFT

    subroutine getFFTWKind(bctype, kind_fwd, kind_bwd)
        implicit none
        character(2), intent(in) :: bctype
        integer(c_int), intent(out) :: kind_fwd, kind_bwd

        select case(bctype)
        case('PP')
            kind_fwd = FFTW_R2HC ! or FFTW_DHT(but slightly slower)
            kind_bwd = FFTW_HC2R ! or FFTW_DHT(but slightly slower)
        case('NN')
            kind_fwd = FFTW_REDFT10
            kind_bwd = FFTW_REDFT01
        case('DD')
            kind_fwd = FFTW_RODFT10
            kind_bwd = FFTW_RODFT01
        case('ND')
            kind_fwd = FFTW_REDFT11
            kind_bwd = FFTW_REDFT11
        case('DN')
            kind_fwd = FFTW_RODFT11
            kind_bwd = FFTW_RODFT11
        end select

        return
    end subroutine getFFTWKind

    subroutine getNormFactor(bctype, normfactor)
        implicit none
        character(2), intent(in) :: bctype
        real(fp), intent(out) :: normfactor

        select case(bctype)
        case('PP')
            normfactor = 1.0_fp
        case('NN')
            normfactor = 2.0_fp
        case('DD')
            normfactor = 2.0_fp
        case('ND')
            normfactor = 2.0_fp
        case('DN')
            normfactor = 2.0_fp
        end select

        return
    end subroutine getNormFactor

    subroutine getEigenvalues(ist, isz, isz_global, bctype, lambda)
        implicit none
        integer,      intent(in) :: ist
        integer,      intent(in) :: isz
        integer,      intent(in) :: isz_global
        character(2), intent(in) :: bctype
        real(fp), dimension(isz), intent(out) :: lambda

        real(fp), allocatable, dimension(:) :: lambda_glb, lambda_aux
        real(fp) :: pi
        integer :: i, ien, n

        pi = acos(-1._fp)
        ien = ist+isz-1

        select case(bctype)
        case('PP')
            n = isz_global
            allocate(lambda_glb(n))
            allocate(lambda_aux(n))
            do i = 1, n
                lambda_aux(i) = 2.0_fp*( cos(2.0_fp*pi*(i-1.0_fp)/n) - 1.0_fp )
            enddo

            select case(mkl_packed_fmt)
            case(MKL_PACKED_FMT_PACK)
                ! e.g. 1, assumes n is even, n = 8 here
                !                       1     2     3     4     5     6     7     8
                ! from FFTW-format   (r[0], r[1], r[2], r[3], r[4], i[3], i[2], i[1])
                !                     e[0], e[1], e[2], e[3], e[4], e[3], e[2], e[1]
                ! to MKL PACK format (r[0], r[1], i[1], r[2], i[2], r[3], i[3], r[4])
                !                     e[0], e[1], e[1], e[2], e[2], e[3], e[3], e[4]
                ! 
                ! e.g. 2, assumes n is odd, n = 7 here
                !                       1     2     3     4     5     6     7
                ! from FFTW-format   (r[0], r[1], r[2], r[3], i[3], i[2], i[1])
                !                     e[0], e[1], e[2], e[3], e[3], e[2], e[1]
                ! to MKL PACK format (r[0], r[1], i[1], r[2], i[2], r[3], i[3])
                !                     e[0], e[1], e[1], e[2], e[2], e[3], e[3]
                lambda_glb(1) = lambda_aux(0 + 1)
                do i = 2, n
                    lambda_glb(i)  = lambda_aux(i/2 + 1)
                enddo

            case(MKL_PACKED_FMT_PERM)
                ! e.g. 1, assumes n is even, n = 8 here
                !                       1     2     3     4     5     6     7     8
                ! from FFTW-format   (r[0], r[1], r[2], r[3], r[4], i[3], i[2], i[1])
                !                     e[0], e[1], e[2], e[3], e[4], e[3], e[2], e[1]
                ! to MKL PERM format (r[0], r[4], r[1], i[1], r[2], i[2], r[3], i[3])
                !                     e[0], e[4], e[1], e[1], e[2], e[2], e[3], e[3]
                ! 
                ! e.g. 2, assumes n is odd, n = 7 here
                !                       1     2     3     4     5     6     7
                ! from FFTW-format   (r[0], r[1], r[2], r[3], i[3], i[2], i[1])
                !                     e[0], e[1], e[2], e[3], e[3], e[2], e[1]
                ! to MKL PERM format (r[0], r[1], i[1], r[2], i[2], r[3], i[3])
                !                     e[0], e[1], e[1], e[2], e[2], e[3], e[3]
                if (mod(n,2) == 0) then
                    lambda_glb(1) = lambda_aux(0 + 1)
                    lambda_glb(2) = lambda_aux(n/2 + 1)
                    do i = 3, n
                        lambda_glb(i) = lambda_aux((i-1)/2 + 1)
                    enddo
                else
                    lambda_glb(1) = lambda_aux(0 + 1)
                    do i = 2, n
                        lambda_glb(i)  = lambda_aux(i/2 + 1)
                    enddo
                endif
                
            case(MKL_PACKED_FMT_PERM2)
                ! e.g. 1, assumes n is even, n = 8 here
                !                       1     2     3     4     5     6     7     8
                ! from FFTW-format   (r[0], r[1], r[2], r[3], r[4], i[3], i[2], i[1])
                !                     e[0], e[1], e[2], e[3], e[4], e[3], e[2], e[1]
                ! to PERM2 format    (r[0], r[4], r[1], i[1], r[2], i[2], r[3], i[3])
                !                     e[0], e[4], e[1], e[1], e[2], e[2], e[3], e[3]
                ! 
                ! e.g. 2, assumes n is odd, n = 7 here
                !                       1     2     3     4     5     6     7
                ! from FFTW-format   (r[0], r[1], r[2], r[3], i[3], i[2], i[1])
                !                     e[0], e[1], e[2], e[3], e[3], e[2], e[1]
                ! to PERM2 format    (r[0], i[3], r[1], i[1], r[2], i[2], r[3])
                !                     e[0], e[3], e[1], e[1], e[2], e[2], e[3]
                lambda_glb(1) = lambda_aux(0 + 1)
                lambda_glb(2) = lambda_aux(n/2 + 1)
                do i = 3, n
                    lambda_glb(i) = lambda_aux((i-1)/2 + 1)
                enddo
                
            end select
            
            do i = ist, ien; lambda(i-ist+1) = lambda_glb(i); enddo
            deallocate(lambda_glb)
            deallocate(lambda_aux)
        case('NN')
            do i = ist, ien; lambda(i-ist+1) = 2.0_fp*( cos(pi*(i-1.0_fp)/isz_global) - 1.0_fp ); enddo
        case('DD')
            do i = ist, ien; lambda(i-ist+1) = 2.0_fp*( cos(pi*(i-0.0_fp)/isz_global) - 1.0_fp ); enddo
        case('ND','DN')
            do i = ist, ien; lambda(i-ist+1) = 2.0_fp*( cos(pi*(2*i-1.0_fp)/(2.*isz_global)) - 1.0_fp ); enddo
        end select

        return
    end subroutine getEigenvalues

    subroutine executeFFT(plan, work)
        implicit none
        type(fft_mkl_plan), intent(inout) :: plan
        real(fp), dimension(*), intent(inout) :: work

        call fft_mkl_execute_r2r(plan, work)

        return
    end subroutine executeFFT

    subroutine freeFFT(fft_plan)
        implicit none
        type(fft_mkl_plan), dimension(2,2), intent(inout) :: fft_plan
        integer :: idir, i
        
        do idir = 1, 2
            do i = 1, 2
                call fft_mkl_destroy_plan(fft_plan(i,idir))
            enddo
        enddo

        return
    end subroutine freeFFT
    
end module mod_fft