module mod_hdf5
    use mod_type, only: fp
    use hdf5
    
    implicit none
    
    include 'mpif.h'

    ! make everything private unless declared public
    private

    integer, parameter :: xfer_size_limit = 2147483647
    integer, parameter :: xfer_size_batch = 1073741824
    ! Each process reads/writes a certain amount of data from/to a subset of a dataset at a time,
    ! if the total size (in byte) of the transfer data is larger than 2^31-1 Bytes (slightly
    ! smaller than 2 GB), then read/write operations will fail because of the limitation of
    ! MPI-IO. Hence, the too large data should be transfered in batches. "nbatch" controls
    ! the number of batches.

    integer, save :: comm, myrank

    integer :: hdferr
    integer :: ks_xfer, ke_xfer
    integer :: rank
    integer(HID_T) :: memspace      ! Dataspace identifier in memory
    integer(HID_T) :: filespace     ! Dataspace identifier in file
    ! integer(HID_T) :: fileid        ! File identifier
    integer(HID_T) :: dsetid        ! Dataset identifier
    integer(HID_T) :: plistid       ! Property list identifier
    integer(HSIZE_T), dimension(1) :: dims_1d
    integer(HSIZE_T), dimension(2) :: dims_2d
    integer(HSIZE_T), dimension(3) :: dims_3d
    integer(HSIZE_T), dimension(3) :: dims_3d_chunk
    integer(HSIZE_T), dimension(3) :: count
    integer(HSSIZE_T),dimension(3) :: offset
    integer(HSIZE_T), dimension(3) :: stride
    integer(HSIZE_T), dimension(3) :: blocksize

    ! public user routines
    public :: HID_T
    public :: initIO, freeIO, createFile, openFile, closeFile, &
              readAttribute, read1d, read3d, &
              writeAttribute, write1d, write3d
    ! public :: read2d, write2d

    interface readAttribute
        module procedure readAttribute_int
        module procedure readAttribute_real
    end interface readAttribute

    interface writeAttribute
        module procedure writeAttribute_int
        module procedure writeAttribute_real
    end interface writeAttribute

    interface write1d
        module procedure write1d_singleproc_0
        module procedure write1d_singleproc
        module procedure write1d_multiproc
    end interface write1d

    ! interface write2d
    !     module procedure write2d_singleproc
    !     module procedure write2d_multiproc
    ! end interface write2d

    interface write3d
        module procedure write3d_col
        module procedure write3d_ind
    end interface write3d

contains
    subroutine initIO(comm_in)
        implicit none
        integer, intent(in) :: comm_in

        integer :: ierr
        
        comm = comm_in
        call MPI_COMM_RANK(comm, myrank, ierr)

        ! Initialize FORTRAN predefined datatypes
        call h5open_f(hdferr)

        return
    end subroutine initIO
    
    subroutine freeIO
        implicit none

        ! Close FORTRAN predefined datatypes
        call h5close_f(hdferr)

        return
    end subroutine freeIO

    subroutine createFile(filename, fileid)
        implicit none
        character(*), intent(in) :: filename
        integer(HID_T), intent(out) :: fileid
        
        ! Setup file access property list with parallel I/O access
        call h5pcreate_f(H5P_FILE_ACCESS_F, plistid, hdferr)
        call h5pset_fapl_mpio_f(plistid, comm, MPI_INFO_NULL, hdferr)

        ! Create the file collectively
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, fileid, hdferr, access_prp = plistid)
        call h5pclose_f(plistid, hdferr)

        return
    end subroutine createFile

    subroutine openFile(filename, fileid)
        implicit none
        character(*), intent(in) :: filename
        integer(HID_T), intent(out) :: fileid

        logical :: alive
        integer :: ierr
        
        inquire(file=filename, exist=alive)
        if (.not. alive) then
            if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.openFile: "//filename//" doesn't exist!"
            call MPI_FINALIZE(ierr)
            stop
        endif
        
        ! Setup file access property list with parallel I/O access
        call h5pcreate_f(H5P_FILE_ACCESS_F, plistid, hdferr)
        call h5pset_fapl_mpio_f(plistid, comm, MPI_INFO_NULL, hdferr)

        ! Open an existing file
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, fileid, hdferr, access_prp = plistid)
        call h5pclose_f(plistid, hdferr)

        return
    end subroutine openFile

    subroutine closeFile(fileid)
        implicit none
        integer(HID_T), intent(in) :: fileid

        ! Close the file
        call h5fclose_f(fileid, hdferr)
  
        return
    end subroutine closeFile

    subroutine readAttribute_int(fileid, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        character(*), intent(in) :: tag
        integer, intent(out) :: var

        ! Create property list for independent dataset read
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
        ! Open an existing dataset
        call h5dopen_f(fileid, tag, dsetid, hdferr)
        ! Read the dataset independently.
        call h5dread_f(dsetid, H5T_NATIVE_INTEGER, var, dims_1d, hdferr, xfer_prp = plistid)
        ! Close the dataset
        call h5dclose_f(dsetid, hdferr)
        ! Close the property list
        call h5pclose_f(plistid, hdferr)

        return
    end subroutine readAttribute_int

    subroutine readAttribute_real(fileid, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        character(*), intent(in) :: tag
        real(fp), intent(out) :: var
        
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
        call h5dopen_f(fileid, tag, dsetid, hdferr)
#ifdef _SINGLE_PREC
        call h5dread_f(dsetid, H5T_NATIVE_REAL, var, dims_1d, hdferr, xfer_prp = plistid)
#else
        call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_1d, hdferr, xfer_prp = plistid)
#endif
        call h5dclose_f(dsetid, hdferr)
        call h5pclose_f(plistid, hdferr)

        return
    end subroutine readAttribute_real

    subroutine writeAttribute_int(fileid, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        character(*), intent(in) :: tag
        integer, intent(in) :: var

        ! create the data space for the dataset
        rank = 1
        dims_1d = 1
        call h5screate_simple_f(rank, dims_1d, filespace, hdferr)
        ! create the dataset with default properties
        call h5dcreate_f(fileid, tag, H5T_NATIVE_INTEGER, filespace, dsetid, hdferr)
        ! Create property list for independent dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
        ! write the dataset
        if (myrank == 0) call h5dwrite_f(dsetid, H5T_NATIVE_INTEGER, var, dims_1d, hdferr, xfer_prp = plistid)
        ! close the property list
        call h5pclose_f(plistid, hdferr)
        ! close the dataset
        call h5dclose_f(dsetid, hdferr)
        ! terminate access to the dataspace
        call h5sclose_f(filespace, hdferr)

        return
    end subroutine writeAttribute_int

    subroutine writeAttribute_real(fileid, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        character(*), intent(in) :: tag
        real(fp), intent(in) :: var

        rank = 1
        dims_1d = 1
        call h5screate_simple_f(rank, dims_1d, filespace, hdferr)
#ifdef _SINGLE_PREC
        call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
#else
        call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
#endif
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
#ifdef _SINGLE_PREC
        if (myrank == 0) call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var, dims_1d, hdferr, xfer_prp = plistid)
#else
        if (myrank == 0) call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_1d, hdferr, xfer_prp = plistid)
#endif
        call h5pclose_f(plistid, hdferr)
        call h5dclose_f(dsetid, hdferr)
        call h5sclose_f(filespace, hdferr)

        return
    end subroutine writeAttribute_real

    subroutine read1d(fileid, is_involved, st, sz, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        logical, intent(in) :: is_involved
        integer, intent(in) :: st, sz
        character(*), intent(in) :: tag
        real(fp), dimension(sz), intent(out) :: var

        call h5dopen_f(fileid, tag, dsetid, hdferr)
        call h5dget_space_f(dsetid, filespace, hdferr)
        rank = 1
        dims_1d = sz
        call h5screate_simple_f(rank, dims_1d, memspace, hdferr)
        offset(1) = st-1
        count(1) = sz
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset(1), count(1), hdferr)
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
        ! call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, hdferr)
#ifdef _SINGLE_PREC
        if (is_involved) call h5dread_f(dsetid, H5T_NATIVE_REAL, var, dims_1d, hdferr, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#else
        if (is_involved) call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_1d, hdferr, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#endif
        call h5pclose_f(plistid, hdferr)
        call h5dclose_f(dsetid, hdferr)
        call h5sclose_f(filespace, hdferr)
        call h5sclose_f(memspace, hdferr)

        return
    end subroutine read1d

    subroutine write1d_singleproc_0(fileid, sz, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        integer, intent(in) :: sz
        character(*), intent(in) :: tag
        real(fp), dimension(sz) :: var

        rank = 1
        dims_1d = sz
        call h5screate_simple_f(rank, dims_1d, filespace, hdferr)
#ifdef _SINGLE_PREC
        call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
#else
        call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
#endif
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
#ifdef _SINGLE_PREC
        if (myrank == 0) call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var, dims_1d, hdferr, xfer_prp = plistid)
#else
        if (myrank == 0) call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_1d, hdferr, xfer_prp = plistid)
#endif
        call h5pclose_f(plistid, hdferr)
        call h5dclose_f(dsetid, hdferr)
        call h5sclose_f(filespace, hdferr)

        return
    end subroutine write1d_singleproc_0

    subroutine write1d_singleproc(fileid, is_involved, sz, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        logical, intent(in) :: is_involved
        integer, intent(in) :: sz
        character(*), intent(in) :: tag
        real(fp), dimension(sz) :: var

        rank = 1
        dims_1d = sz
        call h5screate_simple_f(rank, dims_1d, filespace, hdferr)
#ifdef _SINGLE_PREC
        call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
#else
        call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
#endif
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
#ifdef _SINGLE_PREC
        if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var, dims_1d, hdferr, xfer_prp = plistid)
#else
        if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_1d, hdferr, xfer_prp = plistid)
#endif
        call h5pclose_f(plistid, hdferr)
        call h5dclose_f(dsetid, hdferr)
        call h5sclose_f(filespace, hdferr)

        return
    end subroutine write1d_singleproc

    subroutine write1d_multiproc(fileid, is_involved, sz_global, st, sz, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        logical, intent(in) :: is_involved
        integer, intent(in) :: sz_global, st, sz
        character(*), intent(in) :: tag
        real(fp), dimension(sz) :: var

        rank = 1
        dims_1d = sz_global
        call h5screate_simple_f(rank, dims_1d, filespace, hdferr)
#ifdef _SINGLE_PREC
        call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
#else
        call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
#endif
        dims_1d = sz
        call h5screate_simple_f(rank, dims_1d, memspace, hdferr)
        offset(1) = st-1
        count(1) = sz
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset(1), count(1), hdferr)
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
#ifdef _SINGLE_PREC
        if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var, dims_1d, hdferr, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#else
        if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_1d, hdferr, &
        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#endif
        call h5pclose_f(plistid, hdferr)
        call h5dclose_f(dsetid, hdferr)
        call h5sclose_f(filespace, hdferr)
        call h5sclose_f(memspace, hdferr)

        return
    end subroutine write1d_multiproc

!     subroutine read2d(fileid, is_involved, st, sz, tag, var)
!         implicit none
!         integer(HID_T), intent(in) :: fileid
!         logical, intent(in) :: is_involved
!         integer, dimension(2), intent(in) :: st, sz
!         character(*), intent(in) :: tag
!         real(fp), dimension(sz(1),sz(2)), intent(out) :: var

!         call h5dopen_f(fileid, tag, dsetid, hdferr)
!         call h5dget_space_f(dsetid, filespace, hdferr)
!         rank = 2
!         dims_2d = sz
!         call h5screate_simple_f(rank, dims_2d, memspace, hdferr)
!         offset(1:2) = st-1
!         count(1:2) = sz
!         call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset(1:2), count(1:2), hdferr)
!         call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
!         call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
!         ! call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, hdferr)
! #ifdef _SINGLE_PREC
!         if (is_involved) call h5dread_f(dsetid, H5T_NATIVE_REAL, var, dims_2d, hdferr, &
!         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
! #else
!         if (is_involved) call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_2d, hdferr, &
!         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
! #endif
!         call h5pclose_f(plistid, hdferr)
!         call h5dclose_f(dsetid, hdferr)
!         call h5sclose_f(filespace, hdferr)
!         call h5sclose_f(memspace, hdferr)

!         return
!     end subroutine read2d

!     subroutine write2d_singleproc(fileid, is_involved, sz, tag, var)
!         implicit none
!         integer(HID_T), intent(in) :: fileid
!         logical, intent(in) :: is_involved
!         integer, dimension(2), intent(in) :: sz
!         character(*), intent(in) :: tag
!         real(fp), dimension(sz(1),sz(2)), intent(in) :: var

!         rank = 2
!         dims_2d = sz
!         call h5screate_simple_f(rank, dims_2d, filespace, hdferr)
! #ifdef _SINGLE_PREC
!         call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
! #else
!         call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
! #endif
!         call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
!         call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
! #ifdef _SINGLE_PREC
!         if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var, dims_2d, hdferr, xfer_prp = plistid)
! #else
!         if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_2d, hdferr, xfer_prp = plistid)
! #endif
!         call h5pclose_f(plistid, hdferr)
!         call h5dclose_f(dsetid, hdferr)
!         call h5sclose_f(filespace, hdferr)

!         return
!     end subroutine write2d_singleproc

!     subroutine write2d_multiproc(fileid, is_involved, sz_global, st, sz, tag, var)
!         implicit none
!         integer(HID_T), intent(in) :: fileid
!         logical, intent(in) :: is_involved
!         integer, dimension(2), intent(in) :: sz_global, st, sz
!         character(*), intent(in) :: tag
!         real(fp), dimension(sz(1),sz(2)), intent(in) :: var

!         rank = 2
!         dims_2d = sz_global
!         call h5screate_simple_f(rank, dims_2d, filespace, hdferr)
! #ifdef _SINGLE_PREC
!         call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
! #else
!         call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
! #endif
!         dims_2d = sz
!         call h5screate_simple_f(rank, dims_2d, memspace, hdferr)
!         offset(1:2) = st-1
!         count(1:2) = sz
!         call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset(1:2), count(1:2), hdferr)
!         call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
!         call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
! #ifdef _SINGLE_PREC
!         if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var, dims_2d, hdferr, &
!         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
! #else
!         if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var, dims_2d, hdferr, &
!         file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
! #endif
!         call h5pclose_f(plistid, hdferr)
!         call h5dclose_f(dsetid, hdferr)
!         call h5sclose_f(filespace, hdferr)
!         call h5sclose_f(memspace, hdferr)

!         return
!     end subroutine write2d_multiproc

    subroutine read3d(fileid, st, sz, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        integer, dimension(3), intent(in) :: st, sz
        character(*), intent(in) :: tag
        real(fp), dimension(1:sz(1),1:sz(2),1:sz(3)), intent(out) :: var

        integer :: nbatch = 1
        integer :: ibatch
        integer(8) :: xfer_size

        ! Open an existing dataset
        call h5dopen_f(fileid, tag, dsetid, hdferr)
        
        ! Get the data space for the whole dataset
        call h5dget_space_f(dsetid, filespace, hdferr)
        
        xfer_size = sz(1)*sz(2)*sz(3)*fp
        if (xfer_size > xfer_size_limit) nbatch = ceiling(xfer_size/(1.0_fp*xfer_size_batch))
        dims_3d_chunk = [sz(1), sz(2), ceiling(sz(3)/(1.0_fp*nbatch))]
        rank = 3
        call h5screate_simple_f(rank, dims_3d_chunk, memspace, hdferr)

        ! Create property list for independent/collective dataset read
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        ! call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, hdferr)

        ! Select hyperslab in the file
        stride(1:3) = 1
        count (1:3) = 1
        blocksize(1) = dims_3d_chunk(1)
        blocksize(2) = dims_3d_chunk(2)
        blocksize(3) = dims_3d_chunk(3)
        do ibatch = 1, nbatch
            offset(1) = st(1)-1
            offset(2) = st(2)-1
            offset(3) = st(3)-1 + (ibatch-1)*blocksize(3)
            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, hdferr, stride, blocksize)
            ! Array index range of the 3rd dimension
            ks_xfer = 1 + (ibatch-1)*blocksize(3)
            ke_xfer = ibatch*blocksize(3)
            ! Reset the memspace & filespace size of the last batch
            if (xfer_size > xfer_size_limit .and. ibatch == nbatch) THEN
                ke_xfer = sz(3)
                blocksize(3) = sz(3) - (nbatch-1)*blocksize(3)
                call h5sset_extent_simple_f(memspace, rank, blocksize, blocksize, hdferr)
                call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, hdferr, stride, blocksize)
            endif
            ! Read the dataset 
#ifdef _SINGLE_PREC
            call h5dread_f(dsetid, H5T_NATIVE_REAL, var(:,:,ks_xfer:ke_xfer), blocksize, hdferr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#else
            call h5dread_f(dsetid, H5T_NATIVE_DOUBLE, var(:,:,ks_xfer:ke_xfer), blocksize, hdferr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#endif
        enddo

        ! Close the property list
        call h5pclose_f(plistid, hdferr)
        ! Close the dataset
        call h5dclose_f(dsetid, hdferr)
        ! Close dataspaces
        call h5sclose_f(filespace, hdferr)
        call h5sclose_f(memspace, hdferr)

        return
    end subroutine read3d

    subroutine write3d_col(fileid, sz_global, st, sz, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        integer, dimension(3), intent(in) :: sz_global, st, sz
        character(*), intent(in) :: tag
        real(fp), dimension(1:sz(1),1:sz(2),1:sz(3)), intent(in) :: var

        integer :: nbatch = 1
        integer :: ibatch
        integer(8) :: xfer_size

        ! Create the data space for the whole dataset
        rank = 3
        dims_3d = sz_global
        call h5screate_simple_f(rank, dims_3d, filespace, hdferr)
        
        xfer_size = sz(1)*sz(2)*sz(3)*fp
        if (xfer_size > xfer_size_limit) nbatch = ceiling(xfer_size/(1.0_fp*xfer_size_batch))
        dims_3d_chunk = [sz(1), sz(2), ceiling(sz(3)/(1.0_fp*nbatch))]
        call h5screate_simple_f(rank, dims_3d_chunk, memspace, hdferr)

        ! Create datasets with default properties
#ifdef _SINGLE_PREC
        call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
#else
        call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
#endif

        ! Create property list for independent/collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        ! call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, hdferr)

        ! Select hyperslab in the file
        stride(1:3) = 1
        count (1:3) = 1
        blocksize(1) = dims_3d_chunk(1)
        blocksize(2) = dims_3d_chunk(2)
        blocksize(3) = dims_3d_chunk(3)
        do ibatch = 1, nbatch
            offset(1) = st(1)-1
            offset(2) = st(2)-1
            offset(3) = st(3)-1 + (ibatch-1)*blocksize(3)
            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, hdferr, stride, blocksize)
            ! Array index range of the 3rd dimension
            ks_xfer = 1 + (ibatch-1)*blocksize(3)
            ke_xfer = ibatch*blocksize(3)
            ! Reset the memspace & filespace size of the last batch
            if (xfer_size > xfer_size_limit .and. ibatch == nbatch) THEN
                ke_xfer = sz(3)
                blocksize(3) = sz(3) - (nbatch-1)*blocksize(3)
                call h5sset_extent_simple_f(memspace, rank, blocksize, blocksize, hdferr)
                call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, hdferr, stride, blocksize)
            endif
            ! Write the dataset 
#ifdef _SINGLE_PREC
            call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var(:,:,ks_xfer:ke_xfer), blocksize, hdferr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#else
            call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var(:,:,ks_xfer:ke_xfer), blocksize, hdferr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#endif
        enddo

        ! Close the property list
        call h5pclose_f(plistid, hdferr)
        ! Close the dataset
        call h5dclose_f(dsetid, hdferr)
        ! Close dataspaces
        call h5sclose_f(filespace, hdferr)
        call h5sclose_f(memspace, hdferr)

        return
    end subroutine write3d_col

    subroutine write3d_ind(fileid, is_involved, sz_global, st, sz, tag, var)
        implicit none
        integer(HID_T), intent(in) :: fileid
        logical, intent(in) :: is_involved
        integer, dimension(3), intent(in) :: sz_global, st, sz
        character(*), intent(in) :: tag
        real(fp), dimension(1:sz(1),1:sz(2),1:sz(3)), intent(in) :: var

        integer :: nbatch = 1
        integer :: ibatch
        integer(8) :: xfer_size

        ! Create the data space for the whole dataset
        rank = 3
        dims_3d = sz_global
        call h5screate_simple_f(rank, dims_3d, filespace, hdferr)
        
        xfer_size = sz(1)*sz(2)*sz(3)*fp
        if (xfer_size > xfer_size_limit) nbatch = ceiling(xfer_size/(1.0_fp*xfer_size_batch))
        dims_3d_chunk = [sz(1), sz(2), ceiling(sz(3)/(1.0_fp*nbatch))]
        call h5screate_simple_f(rank, dims_3d_chunk, memspace, hdferr)

        ! Create datasets with default properties
#ifdef _SINGLE_PREC
        call h5dcreate_f(fileid, tag, H5T_NATIVE_REAL, filespace, dsetid, hdferr)
#else
        call h5dcreate_f(fileid, tag, H5T_NATIVE_DOUBLE, filespace, dsetid, hdferr)
#endif

        ! Create property list for independent/collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plistid, hdferr) 
        call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_INDEPENDENT_F, hdferr)
        ! call h5pset_dxpl_mpio_f(plistid, H5FD_MPIO_COLLECTIVE_F, hdferr)

        ! Select hyperslab in the file
        stride(1:3) = 1
        count (1:3) = 1
        blocksize(1) = dims_3d_chunk(1)
        blocksize(2) = dims_3d_chunk(2)
        blocksize(3) = dims_3d_chunk(3)
        do ibatch = 1, nbatch
            offset(1) = st(1)-1
            offset(2) = st(2)-1
            offset(3) = st(3)-1 + (ibatch-1)*blocksize(3)
            call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, hdferr, stride, blocksize)
            ! Array index range of the 3rd dimension
            ks_xfer = 1 + (ibatch-1)*blocksize(3)
            ke_xfer = ibatch*blocksize(3)
            ! Reset the memspace & filespace size of the last batch
            if (xfer_size > xfer_size_limit .and. ibatch == nbatch) THEN
                ke_xfer = sz(3)
                blocksize(3) = sz(3) - (nbatch-1)*blocksize(3)
                call h5sset_extent_simple_f(memspace, rank, blocksize, blocksize, hdferr)
                call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, hdferr, stride, blocksize)
            endif
            ! Write the dataset 
#ifdef _SINGLE_PREC
            if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_REAL, var(:,:,ks_xfer:ke_xfer), blocksize, hdferr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#else
            if (is_involved) call h5dwrite_f(dsetid, H5T_NATIVE_DOUBLE, var(:,:,ks_xfer:ke_xfer), blocksize, hdferr, &
            file_space_id = filespace, mem_space_id = memspace, xfer_prp = plistid)
#endif
        enddo

        ! Close the property list
        call h5pclose_f(plistid, hdferr)
        ! Close the dataset
        call h5dclose_f(dsetid, hdferr)
        ! Close dataspaces
        call h5sclose_f(filespace, hdferr)
        call h5sclose_f(memspace, hdferr)

        return
    end subroutine write3d_ind
    
end module mod_hdf5