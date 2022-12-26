module decomp_2d_c
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), bind(C, name="decomp_2d_nx_global") :: nx_global
    integer(c_int), bind(C, name="decomp_2d_ny_global") :: ny_global
    integer(c_int), bind(C, name="decomp_2d_nz_global") :: nz_global
    integer(c_int), bind(C, name="decomp_2d_nrank") :: nrank
    integer(c_int), bind(C, name="decomp_2d_nproc") :: nproc
    integer(c_int), bind(C, name="decomp_2d_comm_cart_x_f") :: decomp_2d_comm_cart_x
    integer(c_int), bind(C, name="decomp_2d_comm_cart_y_f") :: decomp_2d_comm_cart_y
    integer(c_int), bind(C, name="decomp_2d_comm_cart_z_f") :: decomp_2d_comm_cart_z

    type, bind(C) :: decomp_2d_info
        type(c_ptr) :: x1cnts, x1disp
        type(c_ptr) :: y1cnts, y1disp
        type(c_ptr) :: y2cnts, y2disp
        type(c_ptr) :: z2cnts, z2disp
        integer(c_int) :: x1count, y1count, y2count, z2count
        integer(c_int), dimension(3) :: xst, xen, xsz
        integer(c_int), dimension(3) :: yst, yen, ysz
        integer(c_int), dimension(3) :: zst, zen, zsz
        type(c_ptr) :: x1st, x1en, x1dist
        type(c_ptr) :: y1st, y1en, y1dist
        type(c_ptr) :: y2st, y2en, y2dist
        type(c_ptr) :: z2st, z2en, z2dist
        logical(c_bool) :: even
    end type
    type(decomp_2d_info), bind(C, name="decomp_main") :: decomp_main
    
    interface
        subroutine decomp_2d_init(nx, ny, nz, p_row, p_col, periodic_bc) bind(C, name='decomp_2d_init')
            import
            integer(c_int), value :: nx, ny, nz
            integer(c_int), value :: p_row, p_col
            logical(c_bool), dimension(*), intent(in) :: periodic_bc
        end subroutine decomp_2d_init

        subroutine decomp_2d_finalize() bind(C, name='decomp_2d_finalize')

        end subroutine decomp_2d_finalize

        subroutine decomp_2d_abort(errorcode, msg) bind(C, name='decomp_2d_abort')
            import
            integer(c_int), value :: errorcode
            character(c_char), dimension(*), intent(in) :: msg
        end subroutine decomp_2d_abort

        subroutine decomp_info_init(decomp, nx, ny, nz) bind(C, name='decomp_info_init')
            import
            type(c_ptr), value :: decomp
            integer(c_int) :: nx, ny, nz
        end subroutine decomp_info_init

        subroutine decomp_info_finalize(decomp) bind(C, name='decomp_info_finalize')
            import
            type(c_ptr), value :: decomp
        end subroutine decomp_info_finalize

        subroutine transpose_x_to_y_real(src, dst, src_size, dst_size, decomp) bind(C, name='transpose_x_to_y_real')
            import
            real(c_double), dimension(*), intent(in ) :: src
            real(c_double), dimension(*), intent(out) :: dst
            integer(c_int), dimension(3), intent(in ) :: src_size, dst_size
            type(c_ptr), value :: decomp
        end subroutine transpose_x_to_y_real

        subroutine transpose_y_to_x_real(src, dst, src_size, dst_size, decomp) bind(C, name='transpose_y_to_x_real')
            import
            real(c_double), dimension(*), intent(in ) :: src
            real(c_double), dimension(*), intent(out) :: dst
            integer(c_int), dimension(3), intent(in ) :: src_size, dst_size
            type(c_ptr), value :: decomp
        end subroutine transpose_y_to_x_real

        subroutine transpose_y_to_z_real(src, dst, src_size, dst_size, decomp) bind(C, name='transpose_y_to_z_real')
            import
            real(c_double), dimension(*), intent(in ) :: src
            real(c_double), dimension(*), intent(out) :: dst
            integer(c_int), dimension(3), intent(in ) :: src_size, dst_size
            type(c_ptr), value :: decomp
        end subroutine transpose_y_to_z_real

        subroutine transpose_z_to_y_real(src, dst, src_size, dst_size, decomp) bind(C, name='transpose_z_to_y_real')
            import
            real(c_double), dimension(*), intent(in ) :: src
            real(c_double), dimension(*), intent(out) :: dst
            integer(c_int), dimension(3), intent(in ) :: src_size, dst_size
            type(c_ptr), value :: decomp
        end subroutine transpose_z_to_y_real
    end interface

end module decomp_2d_c