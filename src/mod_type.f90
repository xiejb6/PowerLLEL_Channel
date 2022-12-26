module mod_type

    implicit none

    ! Single precision or double precision (default)
#ifdef _SINGLE_PREC
    integer, parameter :: fp = selected_real_kind(6)
#else
    integer, parameter :: fp = selected_real_kind(15)
#endif

    ! The derived type for recording the info. of statistics processing
    type stat_info_t
        integer :: nts, nte
        integer :: nspl
    end type

end module mod_type