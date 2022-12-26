module mod_calcRHS
    use mod_type,       only: fp
    use mod_parameters, only: dt, nhalo
    use mod_mpi,        only: sz
    use mod_mesh,       only: dx_inv, dy_inv, dzf_inv
    !$ use omp_lib

    implicit none

contains
    subroutine calcRHS(u, v, w, rhs)
        implicit none
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in ) :: u, v, w
        real(fp), dimension(0:,0:,0:), intent(out) :: rhs

        real(fp) :: dtidxi, dtidyi, dti
        integer  :: i, j, k

        dti = 1.0_fp/dt
        dtidxi = dti * dx_inv
        dtidyi = dti * dy_inv

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            rhs(i,j,k) = (u(i,j,k)-u(i-1,j,k))*dtidxi + &
                         (v(i,j,k)-v(i,j-1,k))*dtidyi + &
                         (w(i,j,k)-w(i,j,k-1))*dti*dzf_inv(k)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        return
    end subroutine calcRHS

end module mod_calcRHS