module mod_dataIO
    use mod_type
    use mod_parameters, only: nx, ny, nz, p_row, p_col, nhalo, &
                              nt_out_inst, fn_prefix_input_inst, fn_prefix_inst, &
                              nt_out_save, fn_prefix_save, overwrite_save, auto_cleanup, num_retained, &
                              nt_out_stat, nt_init_stat, fn_prefix_input_stat, fn_prefix_stat, &
                              stat_info, stat_which_var, nt_out_moni
    use mod_mpi,        only: st, sz, ierr, myrank, coord_pen
    use mod_monitor,    only: outputMonitor
    use mod_statistics, only: u_stat, v_stat, w_stat, u2_stat, v2_stat, w2_stat, &
                              uv_stat,uw_stat,vw_stat, p_stat, p2_stat
    use mod_hdf5,       only: HID_T, openFile, createFile, closeFile, readAttribute, writeAttribute, read3d, write3d

    implicit none

    include 'mpif.h'

    ! make everything private unless declared public
    private

    ! public user routines
    public :: inputData, outputData, inputStatData, outputStatData, &
              inputField, outputField, outputFieldWithHalo

contains
    subroutine inputData(nt_in, u, v, w)
        implicit none
        integer, intent(out) :: nt_in
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(out) :: u, v, w

        integer :: nt_last

        if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputData: Reading checkpoint fields ..."

        call inputField(trim(adjustl(fn_prefix_input_inst))//'_u.h5', nt_last, 'u', u(1:sz(1),1:sz(2),1:sz(3)))
        if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputData: Finish reading checkpoint field <u>!"
        nt_in = nt_last
        call inputField(trim(adjustl(fn_prefix_input_inst))//'_v.h5', nt_last, 'v', v(1:sz(1),1:sz(2),1:sz(3)))
        if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputData: Finish reading checkpoint field <v>!"
        call checkTimestamp(nt_in, nt_last, 'v')
        call inputField(trim(adjustl(fn_prefix_input_inst))//'_w.h5', nt_last, 'w', w(1:sz(1),1:sz(2),1:sz(3)))
        if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputData: Finish reading checkpoint field <w>!"
        call checkTimestamp(nt_in, nt_last, 'w')

        return
    contains
        subroutine checkTimestamp(base, actual, vartag)
            integer, intent(in) :: base, actual
            character(*), intent(in) :: vartag
            if (actual /= base) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.inputData: The timestamp of <"//vartag// &
                                                "> does not match with that of <u>!"
                call MPI_FINALIZE(ierr)
                stop
            endif
        end subroutine
    end subroutine inputData

    subroutine outputData(nt, u, v, w, p, u_crf)
        implicit none
        integer, intent(in) :: nt
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in) :: u, v, w
        real(fp), dimension(0:,0:,0:), intent(in) :: p
        real(fp), intent(in) :: u_crf

        character(10) :: string_dump
        integer :: nt_cleanup

        if (mod(nt, nt_out_moni) == 0) then
            call outputMonitor(nt, u, v, w, p, u_crf)
        endif
        
        if (mod(nt, nt_out_inst) == 0) then
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Writing instantaneous fields ..."
            write(string_dump, "('_',I8.8,'_')") nt
            call outputField(fn_prefix_inst//string_dump//'u.h5', nt, 'u', u(1:sz(1),1:sz(2),1:sz(3))+u_crf)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Finish writing inst. field <u>!"
            call outputField(fn_prefix_inst//string_dump//'v.h5', nt, 'v', v(1:sz(1),1:sz(2),1:sz(3)))
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Finish writing inst. field <v>!"
            call outputField(fn_prefix_inst//string_dump//'w.h5', nt, 'w', w(1:sz(1),1:sz(2),1:sz(3)))
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Finish writing inst. field <w>!"
            call outputField(fn_prefix_inst//string_dump//'p.h5', nt, 'p', p(1:sz(1),1:sz(2),1:sz(3)))
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Finish writing inst. field <p>!"
        endif
        
        if (mod(nt, nt_out_save) == 0) then
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Writing checkpoint fields ..."
            if (overwrite_save) then
                string_dump = '_'
            else
                write(string_dump, "('_',I8.8,'_')") nt
            endif
            call outputField(fn_prefix_save//trim(string_dump)//'u.h5', nt, 'u', u(1:sz(1),1:sz(2),1:sz(3))+u_crf)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Finish writing checkpoint field <u>!"
            call outputField(fn_prefix_save//trim(string_dump)//'v.h5', nt, 'v', v(1:sz(1),1:sz(2),1:sz(3)))
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Finish writing checkpoint field <v>!"
            call outputField(fn_prefix_save//trim(string_dump)//'w.h5', nt, 'w', w(1:sz(1),1:sz(2),1:sz(3)))
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputData: Finish writing checkpoint field <w>!"
        endif

        if (nt > nt_init_stat) then
            call outputStatData(nt)
        endif

        if (myrank == 0) then
        if ((.not.overwrite_save) .and. auto_cleanup .and. mod(nt, nt_out_save) == 0) then
            nt_cleanup = nt-(num_retained+1)*nt_out_save
            write(string_dump, "('_',I8.8,'_')") nt_cleanup
            if (nt_cleanup > 0) then
                write(*,'(A)') "PowerLLEL.NOTE.outputData: Automatically cleanup of inst. checkpoint files ..."
                call system("rm "//fn_prefix_save//trim(string_dump)//'u.h5')
                call system("rm "//fn_prefix_save//trim(string_dump)//'v.h5')
                call system("rm "//fn_prefix_save//trim(string_dump)//'w.h5')
            endif
            if (nt_cleanup > nt_init_stat) then
                write(*,'(A)') "PowerLLEL.NOTE.outputData: Automatically cleanup of stat. checkpoint files ..."
                if (stat_which_var(1)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_u.h5')
                if (stat_which_var(2)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_v.h5')
                if (stat_which_var(3)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_w.h5')
                if (stat_which_var(4)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_u2.h5')
                if (stat_which_var(5)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_v2.h5')
                if (stat_which_var(6)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_w2.h5')
                if (stat_which_var(7)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_uv.h5')
                if (stat_which_var(8)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_uw.h5')
                if (stat_which_var(9)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_vw.h5')
                if (stat_which_var(10)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_p.h5')
                if (stat_which_var(11)) &
                    call system("rm "//fn_prefix_save//trim(string_dump)//'stat_p2.h5')
            endif
        endif
        endif

        ! call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        return
    end subroutine outputData

    subroutine inputStatData(nt_in_inst)
        implicit none
        integer, intent(in) :: nt_in_inst

        integer :: nt_last

        if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Reading checkpoint fields of statistics ..."

        if (stat_which_var( 1)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_u.h5', nt_last, 'u_stat', u_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <u_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'u_stat')
        endif
        if (stat_which_var( 2)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_v.h5', nt_last, 'v_stat', v_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <v_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'v_stat')
        endif
        if (stat_which_var( 3)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_w.h5', nt_last, 'w_stat', w_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <w_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'w_stat')
        endif
        if (stat_which_var( 4)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_u2.h5', nt_last, 'u2_stat', u2_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <u2_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'u2_stat')
        endif
        if (stat_which_var( 5)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_v2.h5', nt_last, 'v2_stat', v2_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <v2_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'v2_stat')
        endif
        if (stat_which_var( 6)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_w2.h5', nt_last, 'w2_stat', w2_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <w2_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'w2_stat')
        endif
        if (stat_which_var( 7)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_uv.h5', nt_last, 'uv_stat', uv_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <uv_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'uv_stat')
        endif
        if (stat_which_var( 8)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_uw.h5', nt_last, 'uw_stat', uw_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <uw_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'uw_stat')
        endif
        if (stat_which_var( 9)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_vw.h5', nt_last, 'vw_stat', vw_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <vw_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'vw_stat')
        endif
        if (stat_which_var(10)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_p.h5', nt_last, 'p_stat', p_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <p_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'p_stat')
        endif
        if (stat_which_var(11)) then
            call inputField(trim(adjustl(fn_prefix_input_stat))//'_p2.h5', nt_last, 'p2_stat', p2_stat, stat_info)
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.inputStatData: Finish reading checkpoint field <p2_stat>!"
            call checkTimestamp(nt_in_inst, nt_last, 'p2_stat')
        endif

        return
    contains
        subroutine checkTimestamp(base, actual, vartag)
            integer, intent(in) :: base, actual
            character(*), intent(in) :: vartag
            if (actual /= base) then
                if (myrank == 0) write(*,'(A)') "PowerLLEL.ERROR.inputStatData: The timestamp of <"//vartag// &
                                                "> does not match with that of inst. field!"
                call MPI_FINALIZE(ierr)
                stop
            endif
        end subroutine
    end subroutine inputStatData

    subroutine outputStatData(nt)
        implicit none
        integer, intent(in) :: nt

        character(10) :: string_dump
        character(19) :: string_dump2

        if (mod(nt, nt_out_save) == 0) then
            
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Writing checkpoint fields of statistics ..."
            if (overwrite_save) then
                string_dump = '_'
            else
                write(string_dump, "('_',I8.8,'_')") nt
            endif
            if (stat_which_var( 1)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_u.h5', nt, 'u_stat',  u_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <u_stat>!"
            endif
            if (stat_which_var( 2)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_v.h5', nt, 'v_stat',  v_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <v_stat>!"
            endif
            if (stat_which_var( 3)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_w.h5', nt, 'w_stat',  w_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <w_stat>!"
            endif
            if (stat_which_var( 4)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_u2.h5', nt, 'u2_stat', u2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <u2_stat>!"
            endif
            if (stat_which_var( 5)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_v2.h5', nt, 'v2_stat', v2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <v2_stat>!"
            endif
            if (stat_which_var( 6)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_w2.h5', nt, 'w2_stat', w2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <w2_stat>!"
            endif
            if (stat_which_var( 7)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_uv.h5', nt, 'uv_stat', uv_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <uv_stat>!"
            endif
            if (stat_which_var( 8)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_uw.h5', nt, 'uw_stat', uw_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <uw_stat>!"
            endif
            if (stat_which_var( 9)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_vw.h5', nt, 'vw_stat', vw_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <vw_stat>!"
            endif
            if (stat_which_var(10)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_p.h5', nt, 'p_stat',  p_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <p_stat>!"
            endif
            if (stat_which_var(11)) then
                call outputField(fn_prefix_save//trim(string_dump)//'stat_p2.h5', nt, 'p2_stat', p2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing checkpoint field <p2_stat>!"
            endif

        endif

        if (mod(nt-nt_init_stat,nt_out_stat) == 0) then
            
            if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Writing statistics fields ..."
            
            write(string_dump2,"('_',I8.8,'-',I8.8,'_')") stat_info%nts, stat_info%nte
            if (stat_which_var( 1)) then
                call outputField(fn_prefix_stat//string_dump2//'u.h5', nt, 'u_stat',  u_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <u_stat>!"
            endif
            if (stat_which_var( 2)) then
                call outputField(fn_prefix_stat//string_dump2//'v.h5', nt, 'v_stat',  v_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <v_stat>!"
            endif
            if (stat_which_var( 3)) then
                call outputField(fn_prefix_stat//string_dump2//'w.h5', nt, 'w_stat',  w_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <w_stat>!"
            endif
            if (stat_which_var( 4)) then
                call outputField(fn_prefix_stat//string_dump2//'u2.h5', nt, 'u2_stat', u2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <u2_stat>!"
            endif
            if (stat_which_var( 5)) then
                call outputField(fn_prefix_stat//string_dump2//'v2.h5', nt, 'v2_stat', v2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <v2_stat>!"
            endif
            if (stat_which_var( 6)) then
                call outputField(fn_prefix_stat//string_dump2//'w2.h5', nt, 'w2_stat', w2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <w2_stat>!"
            endif
            if (stat_which_var( 7)) then
                call outputField(fn_prefix_stat//string_dump2//'uv.h5', nt, 'uv_stat', uv_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <uv_stat>!"
            endif
            if (stat_which_var( 8)) then
                call outputField(fn_prefix_stat//string_dump2//'uw.h5', nt, 'uw_stat', uw_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <uw_stat>!"
            endif
            if (stat_which_var( 9)) then
                call outputField(fn_prefix_stat//string_dump2//'vw.h5', nt, 'vw_stat', vw_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <vw_stat>!"
            endif
            if (stat_which_var(10)) then
                call outputField(fn_prefix_stat//string_dump2//'p.h5', nt, 'p_stat',  p_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <p_stat>!"
            endif
            if (stat_which_var(11)) then
                call outputField(fn_prefix_stat//string_dump2//'p2.h5', nt, 'p2_stat', p2_stat, stat_info)
                if (myrank == 0) write(*,'(A)') "PowerLLEL.NOTE.outputStatData: Finish writing stat. field <p2_stat>!"
            endif

        endif

        return
    end subroutine outputStatData

    subroutine inputField(fn_field, nt, vartag, var, stat_info)
        implicit none
        character(*), intent(in) :: fn_field
        integer, intent(out) :: nt
        character(*), intent(in) :: vartag
        real(fp), dimension(1:sz(1),1:sz(2),1:sz(3)), intent(out) :: var
        type(stat_info_t), intent(out), optional :: stat_info

        integer(HID_T) :: fh

        if (.not. present(stat_info)) then
            call openFile(fn_field, fh)
            call readAttribute(fh, 'nt', nt)
            call read3d(fh, st, sz, vartag, var)
            call closeFile(fh)
        else
            call openFile(fn_field, fh)
            call readAttribute(fh, 'nts',  stat_info%nts )
            call readAttribute(fh, 'nte',  stat_info%nte )
            call readAttribute(fh, 'nspl', stat_info%nspl)
            call read3d(fh, st, sz, vartag, var)
            call closeFile(fh)

            nt = stat_info%nte
        endif

        return
    end subroutine inputField

    subroutine outputField(fn_field, nt, vartag, var, stat_info)
        implicit none
        character(*), intent(in) :: fn_field
        integer, intent(in) :: nt
        character(*), intent(in) :: vartag
        real(fp), dimension(1:sz(1),1:sz(2),1:sz(3)), intent(in) :: var
        type(stat_info_t), intent(in), optional :: stat_info

        integer(HID_T) :: fh

        if (.not. present(stat_info)) then
            call createFile(fn_field, fh)
            call writeAttribute(fh, 'nt', nt)
            call write3d(fh, (/nx,ny,nz/), st, sz, vartag, var)
            call closeFile(fh)
        else
            call createFile(fn_field, fh)
            call writeAttribute(fh, 'nts',  stat_info%nts )
            call writeAttribute(fh, 'nte',  stat_info%nte )
            call writeAttribute(fh, 'nspl', stat_info%nspl)
            call write3d(fh, (/nx,ny,nz/), st, sz, vartag, var)
            call closeFile(fh)
        end if

        return
    end subroutine outputField

    subroutine outputFieldWithHalo(fn_field, nt, nly, vartag, var)
        implicit none
        character(*), intent(in) :: fn_field
        integer, intent(in) :: nt
        integer, dimension(6), intent(in) :: nly
        character(*), intent(in) :: vartag
        real(fp), dimension(1-nly(1):sz(1)+nly(2),1-nly(3):sz(2)+nly(4),1-nly(5):sz(3)+nly(6)), intent(in) :: var

        integer(HID_T) :: fh
        integer :: nx_h, ny_h, nz_h
        integer, dimension(3) :: st_h, sz_h

        nx_h = nx+(nly(1)+nly(2))
        ny_h = ny+(nly(3)+nly(4))*p_row
        nz_h = nz+(nly(5)+nly(6))*p_col
        sz_h(1) = sz(1) + nly(1)+nly(2)
        sz_h(2) = sz(2) + nly(3)+nly(4)
        sz_h(3) = sz(3) + nly(5)+nly(6)
        st_h(1) = 1
        st_h(2) = 1+coord_pen(1)*sz_h(2)
        st_h(3) = 1+coord_pen(2)*sz_h(3)

        call createFile(fn_field, fh)
        call writeAttribute(fh, 'nt', nt)
        call write3d(fh, (/nx_h,ny_h,nz_h/), st_h, sz_h, vartag, var)
        call closeFile(fh)

        return
    end subroutine outputFieldWithHalo

end module mod_dataIO