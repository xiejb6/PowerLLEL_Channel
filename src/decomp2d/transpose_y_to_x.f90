!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contains the routines that transpose data from Y to X pencil

  subroutine transpose_y_to_x_real(src, dst, opt_decomp)

    implicit none
    
    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

#ifdef GPTL
    ret = gptlstart('----mem_split_yx')
#endif

    ! rearrange source array as send buffer
    call mem_split_yx_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%y1dist, decomp)

#ifdef GPTL
    ret = gptlstop('----mem_split_yx')
#endif

#ifdef GPTL
    ret = gptlstart('----mpi_alltoall_yx')
#endif
    
    ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%y1count, &
         real_type, work2_r, decomp%x1count, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#endif

#ifdef GPTL
    ret = gptlstop('----mpi_alltoall_yx')
    ret = gptlstart('----mem_merge_yx')
#endif

    ! rearrange receive buffer
    call mem_merge_yx_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)

#ifdef GPTL
    ret = gptlstop('----mem_merge_yx')
#endif
    
    return
  end subroutine transpose_y_to_x_real

#if 0
  subroutine transpose_y_to_x_complex(src, dst, opt_decomp)

    implicit none
    
    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    
    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)
    
    ! rearrange source array as send buffer
    call mem_split_yx_complex(src, s1, s2, s3, work1_c, dims(1), &
         decomp%y1dist, decomp)
    
    ! transpose using MPI_ALLTOALL(V)
#ifdef EVEN
    call MPI_ALLTOALL(work1_c, decomp%y1count, &
         complex_type, work2_c, decomp%x1count, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_c, decomp%y1cnts, decomp%y1disp, &
         complex_type, work2_c, decomp%x1cnts, decomp%x1disp, &
         complex_type, DECOMP_2D_COMM_COL, ierror)
#endif

    ! rearrange receive buffer
    call mem_merge_yx_complex(work2_c, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)

    return
  end subroutine transpose_y_to_x_complex
#endif

  ! pack/unpack ALLTOALL(V) buffers

  subroutine mem_split_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos
    integer :: k1,k2
    integer :: offset = 0
    !$ integer :: tid, nt, ksz

    k1 = 1
    k2 = n3
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,k,j,i,i1,i2,pos,tid,nt,k1,k2,ksz,offset)
    !$ tid = OMP_GET_THREAD_NUM()
    !$ nt  = OMP_GET_NUM_THREADS()
    !$ call distribute_among_threads('f', nt, tid, n3, 1, ksz, k1, k2)
    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

       !$ offset = (i2-i1+1) * n1 * (k1-1)
#ifdef EVEN
       pos = m * decomp%y1count + 1 + offset
#else
       pos = decomp%y1disp(m) + 1 + offset
#endif

       do k=k1,k2
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL

    return
  end subroutine mem_split_yx_real

#if 0
  subroutine mem_split_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(n1,n2,n3), intent(IN) :: in
    complex(mytype), dimension(*), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2,pos

    do m=0,iproc-1
       if (m==0) then 
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef EVEN
       pos = m * decomp%y1count + 1
#else
       pos = decomp%y1disp(m) + 1
#endif

       do k=1,n3
          do j=i1,i2
             do i=1,n1
                out(pos) = in(i,j,k)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_split_yx_complex
#endif

  subroutine mem_merge_yx_real(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    real(mytype), dimension(*), intent(IN) :: in
    real(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos
    integer :: k1,k2
    integer :: offset = 0
    !$ integer :: tid, nt, ksz

    k1 = 1
    k2 = n3

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,k,j,i,i1,i2,pos,tid,nt,k1,k2,ksz,offset)
    !$ tid = OMP_GET_THREAD_NUM()
    !$ nt  = OMP_GET_NUM_THREADS()
    !$ call distribute_among_threads('f', nt, tid, n3, 1, ksz, k1, k2)
    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

       !$ offset = (i2-i1+1) * n2 * (k1-1)
#ifdef EVEN
       pos = m * decomp%x1count + 1 + offset
#else
       pos = decomp%x1disp(m) + 1 + offset
#endif

       do k=k1,k2
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL

    return
  end subroutine mem_merge_yx_real

#if 0
  subroutine mem_merge_yx_complex(in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none
    
    integer, intent(IN) :: n1,n2,n3
    complex(mytype), dimension(*), intent(IN) :: in
    complex(mytype), dimension(n1,n2,n3), intent(OUT) :: out
    integer, intent(IN) :: iproc
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    
    integer :: i,j,k, m,i1,i2, pos

    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       end if

#ifdef EVEN
       pos = m * decomp%x1count + 1
#else
       pos = decomp%x1disp(m) + 1
#endif

       do k=1,n3
          do j=1,n2
             do i=i1,i2
                out(i,j,k) = in(pos)
                pos = pos + 1
             end do
          end do
       end do
    end do

    return
  end subroutine mem_merge_yx_complex
#endif
