module mod_mpi
#ifdef USE_C
    use, intrinsic :: iso_c_binding, only: c_bool
    use decomp_2d_c, only: decomp_2d_init, decomp_2d_finalize, decomp_main
    use decomp_2d_c, only: decomp_2d_comm_cart_x, decomp_2d_comm_cart_y, decomp_2d_comm_cart_z
#else
    use decomp_2d
#endif
    
    implicit none

    include 'mpif.h'

    ! make everything public unless declared private
    public

#ifdef _SINGLE_PREC
    integer, parameter :: REAL_FP = kind(1.0)
    integer, parameter :: MPI_REAL_FP = MPI_REAL
#else
    integer, parameter :: REAL_FP = kind(1.0d0)
    integer, parameter :: MPI_REAL_FP = MPI_DOUBLE_PRECISION
#endif

    integer, save :: myrank, ierr
    integer, save :: comm_cart
    integer, save :: comm_cart_xpen, comm_cart_ypen, comm_cart_zpen
    integer, dimension(6),   save :: halotype_vel, halotype_one
#ifdef NB_HALO
    ! 1: south
    ! 2: north
    ! 3: bottom
    ! 4: top
    ! 5: south_bottom
    ! 6: south_top
    ! 7: north_bottom
    ! 8: north_top
    integer, dimension(8),   save :: mpitype_nbhalo_vel
    integer, dimension(8),   save :: neighbor_nbhalo
#endif
    integer, dimension(6,3), save :: neighbor_xyz
    integer, dimension(6),   save :: neighbor
    integer, dimension(2),   save :: coord_xpen, coord_ypen, coord_zpen
    integer, dimension(2),   save :: coord_pen
    integer, dimension(3),   save ::  st,  en,  sz
    integer, dimension(3),   save :: xst, xen, xsz
    integer, dimension(3),   save :: yst, yen, ysz
    integer, dimension(3),   save :: zst, zen, zsz

    contains
    subroutine initMPI(nx, ny, nz, bctype_p, p_row, p_col, nhalo)
        implicit none
        integer, intent(in) :: nx, ny, nz
        integer, intent(in) :: p_row, p_col
        character(2), dimension(3), intent(in) :: bctype_p
        integer, dimension(6), intent(in) :: nhalo

        logical, dimension(3) :: periodic_bc
        integer :: bcount, bsize, bstride, oldtype

        periodic_bc(1:3) = .false.
        if ( bctype_p(1) == 'PP' ) periodic_bc(1) = .true.
        if ( bctype_p(2) == 'PP' ) periodic_bc(2) = .true.
        if ( bctype_p(3) == 'PP' ) periodic_bc(3) = .true.
        
#ifdef USE_C
        call decomp_2d_init(nx, ny, nz, p_row, p_col, logical(periodic_bc, c_bool))
#else
        call decomp_2d_init(nx, ny, nz, p_row, p_col, periodic_bc)
#endif
        
        ! staring/ending index and size of data held by current processor       
#ifdef USE_C
        xst = decomp_main%xst; xen = decomp_main%xen; xsz = decomp_main%xsz
        yst = decomp_main%yst; yen = decomp_main%yen; ysz = decomp_main%ysz
        zst = decomp_main%zst; zen = decomp_main%zen; zsz = decomp_main%zsz
#else
        xst(:) = xstart(:); xen(:) = xend(:); xsz(:) = xsize(:)  ! x-pencil
        yst(:) = ystart(:); yen(:) = yend(:); ysz(:) = ysize(:)  ! y-pencil
        zst(:) = zstart(:); zen(:) = zend(:); zsz(:) = zsize(:)  ! z-pencil
#endif

        comm_cart_xpen = DECOMP_2D_COMM_CART_X
        comm_cart_ypen = DECOMP_2D_COMM_CART_Y
        comm_cart_zpen = DECOMP_2D_COMM_CART_Z
        
        ! Find the MPI ranks of neighboring pencils
        !  first dimension 1=west, 2=east, 3=south, 4=north, 5=bottom, 6=top
        ! second dimension 1=x-pencil, 2=y-pencil, 3=z-pencil
        ! x-pencil
        neighbor_xyz(1,1) = MPI_PROC_NULL
        neighbor_xyz(2,1) = MPI_PROC_NULL
        call MPI_CART_SHIFT(comm_cart_xpen, 0, 1, neighbor_xyz(3,1), neighbor_xyz(4,1), ierr)
        call MPI_CART_SHIFT(comm_cart_xpen, 1, 1, neighbor_xyz(5,1), neighbor_xyz(6,1), ierr)
        ! y-pencil
        call MPI_CART_SHIFT(comm_cart_ypen, 0, 1, neighbor_xyz(1,2), neighbor_xyz(2,2), ierr)
        neighbor_xyz(3,2) = MPI_PROC_NULL
        neighbor_xyz(4,2) = MPI_PROC_NULL
        call MPI_CART_SHIFT(comm_cart_ypen, 1, 1, neighbor_xyz(5,2), neighbor_xyz(6,2), ierr)
        ! z-pencil
        call MPI_CART_SHIFT(comm_cart_zpen, 0, 1, neighbor_xyz(1,3), neighbor_xyz(2,3), ierr)
        call MPI_CART_SHIFT(comm_cart_zpen, 1, 1, neighbor_xyz(3,3), neighbor_xyz(4,3), ierr)
        neighbor_xyz(5,3) = MPI_PROC_NULL
        neighbor_xyz(6,3) = MPI_PROC_NULL
        
        ! the coordinary of each process in the x-pencil decomposition
        call MPI_CART_COORDS(comm_cart_xpen, myrank, 2, coord_xpen, ierr)
        coord_ypen = coord_xpen
        coord_zpen = coord_xpen

        ! x-pencil
        st(:) = xst(:); en(:) = xen(:); sz(:) = xsz(:)
        comm_cart = comm_cart_xpen
        neighbor(:) = neighbor_xyz(:,1)
        coord_pen(:) = coord_xpen(:)

#ifdef NB_HALO
        call getNeighborRank2DCart(comm_cart, neighbor_nbhalo)
        call createNBHaloMPIType(nhalo, sz, MPI_REAL_FP, mpitype_nbhalo_vel)
        ! write(*,'(A,I3,A,2I3,A,8I3)') '>>> myrank = ', myrank, ', coord = ', coord_pen, ', nb_neighbor = ', neighbor_nbhalo
#endif
        call createHaloMPIType(nhalo, sz, MPI_REAL_FP, halotype_vel)
        call createHaloMPIType((/1,1,1,1,1,1/), sz, MPI_REAL_FP, halotype_one)

        return
    end subroutine initMPI

    subroutine createHaloMPIType(nhalo, sz, oldtype, halotype)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        integer, intent(in) :: oldtype
        integer, dimension(6), intent(out) :: halotype

        integer :: bcount, bsize, bstride, ierr

        ! halo comm in the west/east direction
        bcount  = (sz(2)+nhalo(3)+nhalo(4)) * (sz(3)+nhalo(5)+nhalo(6))
        bsize   = nhalo(1)
        bstride = sz(1)+nhalo(1)+nhalo(2)
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(1), ierr)
        call MPI_TYPE_COMMIT(halotype(1), ierr)
        bsize   = nhalo(2)
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(2), ierr)
        call MPI_TYPE_COMMIT(halotype(2), ierr)
        
        ! halo comm in the south/north direction
        bcount  = sz(3)+nhalo(5)+nhalo(6)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(3)
        bstride = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4))
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(3), ierr)
        call MPI_TYPE_COMMIT(halotype(3), ierr)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(4)
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(4), ierr)
        call MPI_TYPE_COMMIT(halotype(4), ierr)
        
        ! halo comm in the bottom/top direction
        bcount  = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4)) * nhalo(5)
        bsize   = 1
        bstride = 1
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(5), ierr)
        call MPI_TYPE_COMMIT(halotype(5), ierr)
        bcount  = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4)) * nhalo(6)
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(6), ierr)
        call MPI_TYPE_COMMIT(halotype(6), ierr)

        return
    end subroutine createHaloMPIType

    subroutine getNeighborRank2DCart(comm, neighbor)
        implicit none
        integer, intent(in) :: comm
        integer, dimension(8), intent(out) :: neighbor
        
        integer, dimension(2) :: dims, coords, coords_nb
        integer, dimension(-1:1) :: xcoords_nb, ycoords_nb
        logical, dimension(2) :: periods
        integer :: i, j, ierr

        call MPI_CART_GET(comm, 2, dims, periods, coords, ierr)
        
        xcoords_nb(-1) = coords(1)-1
        xcoords_nb( 0) = coords(1)
        xcoords_nb( 1) = coords(1)+1
        if (periods(1)) then
            xcoords_nb(-1) = modulo(xcoords_nb(-1), dims(1))
            xcoords_nb( 1) = modulo(xcoords_nb( 1), dims(1))
        endif
        
        ycoords_nb(-1) = coords(2)-1
        ycoords_nb( 0) = coords(2)
        ycoords_nb( 1) = coords(2)+1
        if (periods(2)) then
            ycoords_nb(-1) = modulo(ycoords_nb(-1), dims(2))
            ycoords_nb( 1) = modulo(ycoords_nb( 1), dims(2))
        endif

        coords_nb(1) = xcoords_nb(-1)
        coords_nb(2) = ycoords_nb( 0)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(1))
        coords_nb(1) = xcoords_nb( 1)
        coords_nb(2) = ycoords_nb( 0)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(2))
        coords_nb(1) = xcoords_nb( 0)
        coords_nb(2) = ycoords_nb(-1)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(3))
        coords_nb(1) = xcoords_nb( 0)
        coords_nb(2) = ycoords_nb( 1)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(4))
        coords_nb(1) = xcoords_nb(-1)
        coords_nb(2) = ycoords_nb(-1)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(5))
        coords_nb(1) = xcoords_nb(-1)
        coords_nb(2) = ycoords_nb( 1)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(6))
        coords_nb(1) = xcoords_nb( 1)
        coords_nb(2) = ycoords_nb(-1)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(7))
        coords_nb(1) = xcoords_nb( 1)
        coords_nb(2) = ycoords_nb( 1)
        call getMpiRankFromCoords2DCart(comm, coords_nb, dims, neighbor(8))

        return
    end subroutine getNeighborRank2DCart

    subroutine getMpiRankFromCoords2DCart(comm, coords, dims, rank)
        implicit none
        integer, intent(in) :: comm
        integer, dimension(2), intent(in) :: coords, dims
        integer, intent(out) :: rank
        
        integer :: ierr

        if( coords(1)<0 .or. coords(1)>=dims(1) .or. coords(2)<0 .or. coords(2)>=dims(2)) then
            rank = MPI_PROC_NULL
        else
            call MPI_CART_RANK(comm, coords, rank, ierr)
        endif

        return
    end subroutine getMpiRankFromCoords2DCart

    subroutine createNBHaloMPIType(nhalo, sz, oldtype, halotype)
        implicit none
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        integer, intent(in) :: oldtype
        integer, dimension(8), intent(out) :: halotype

        integer :: bcount, bsize, bstride, ierr
        
        ! halo exchange in the south/north direction
        bcount  = sz(3)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(3)
        bstride = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4))
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(1), ierr)
        call MPI_TYPE_COMMIT(halotype(1), ierr)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(4)
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(2), ierr)
        call MPI_TYPE_COMMIT(halotype(2), ierr)
        
        ! halo exchange in the bottom/top direction
        bcount  = 1
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * sz(2) * nhalo(5)
        bstride = 1
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(3), ierr)
        call MPI_TYPE_COMMIT(halotype(3), ierr)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * sz(2) * nhalo(6)
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(4), ierr)
        call MPI_TYPE_COMMIT(halotype(4), ierr)

        ! halo exchange in the south_bottom direction
        bcount  = nhalo(5)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(3)
        bstride = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4))
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(5), ierr)
        call MPI_TYPE_COMMIT(halotype(5), ierr)

        ! halo exchange in the south_top direction
        bcount  = nhalo(6)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(3)
        bstride = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4))
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(6), ierr)
        call MPI_TYPE_COMMIT(halotype(6), ierr)

        ! halo exchange in the north_bottom direction
        bcount  = nhalo(5)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(4)
        bstride = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4))
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(7), ierr)
        call MPI_TYPE_COMMIT(halotype(7), ierr)

        ! halo exchange in the north_top direction
        bcount  = nhalo(6)
        bsize   = (sz(1)+nhalo(1)+nhalo(2)) * nhalo(4)
        bstride = (sz(1)+nhalo(1)+nhalo(2)) * (sz(2)+nhalo(3)+nhalo(4))
        call MPI_TYPE_VECTOR(bcount, bsize, bstride, oldtype, halotype(8), ierr)
        call MPI_TYPE_COMMIT(halotype(8), ierr)

        return
    end subroutine createNBHaloMPIType

    subroutine freeMPI()
        implicit none
        integer :: i

        call decomp_2d_finalize()

        do i = 1, 6
            call MPI_TYPE_FREE(halotype_vel(i), ierr)
            call MPI_TYPE_FREE(halotype_one(i), ierr)
        enddo
#ifdef NB_HALO
        do i = 1, 8
            call MPI_TYPE_FREE(mpitype_nbhalo_vel(i), ierr)
        enddo
#endif

        return
    end subroutine freeMPI

end module mod_mpi