module mod_vortex_postproc
    use mod_type,            only: fp
    use mod_param_postproc,  only: nx, ny, nz, nhalo_zero
    use mod_variables_postproc
    use mod_mesh_postproc,   only: dx_inv, dy_inv, dzc_inv, dzf_inv
    use mod_dataio_postproc, only: outputField
    use mod_mpi,             only: st, sz, neighbor, myrank
    !$ use omp_lib
    
    implicit none

    include 'mpif.h'
    
    ! make everything private unless declared public
    private

    public :: calcVelGradTensor, outVorticity, outQ, outLambda2

contains

    subroutine calcVelGradTensor()
        implicit none
        real(fp) :: weightp, weightm
        integer :: i, j, k

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k, weightp, weightm)
        do k = 1, sz(3)
            if (neighbor(5) == MPI_PROC_NULL .and. k == 1) then
                weightp = 0.5*dzc_inv(k)
                weightm = 1.0*dzc_inv(k-1)
            else if (neighbor(6) == MPI_PROC_NULL .and. k == sz(3)) then
                weightp = 1.0*dzc_inv(k)
                weightm = 0.5*dzc_inv(k-1)
            else
                weightp = 0.5*dzc_inv(k)
                weightm = 0.5*dzc_inv(k-1)
            endif
        do j = 1, sz(2)
        do i = 1, sz(1)
            vel_grad%ux(i, j, k) = (u(i, j, k) - u(i-1, j, k))*dx_inv
            vel_grad%uy(i, j, k) = (u(i, j+1, k) - u(i, j-1, k) + u(i-1, j+1, k) - u(i-1, j-1, k))*0.25*dy_inv
            vel_grad%uz(i, j, k) = (u(i, j, k+1) + u(i-1, j, k+1) - u(i, j, k) - u(i-1, j, k))*0.5*weightp + &
                                   (u(i, j, k) + u(i-1, j, k) - u(i, j, k-1) - u(i-1, j, k-1))*0.5*weightm
            vel_grad%vx(i, j, k) = (v(i+1, j, k) - v(i-1, j, k) + v(i+1, j-1, k) - v(i-1, j-1, k))*0.25*dx_inv
            vel_grad%vy(i, j, k) = (v(i, j, k) - v(i, j-1, k))*dy_inv
            vel_grad%vz(i, j, k) = (v(i, j, k+1) + v(i, j-1, k+1) - v(i, j, k) - v(i, j-1, k))*0.5*weightp + &
                                   (v(i, j, k) + v(i, j-1, k) - v(i, j, k-1) - v(i, j-1, k-1))*0.5*weightm
            vel_grad%wx(i, j, k) = (w(i+1, j, k) - w(i-1, j, k) + w(i+1, j, k-1) - w(i-1, j, k-1))*0.25*dx_inv
            vel_grad%wy(i, j, k) = (w(i, j+1, k) - w(i, j-1, k) + w(i, j+1, k-1) - w(i, j-1, k-1))*0.25*dy_inv
            vel_grad%wz(i, j, k) = (w(i, j, k) - w(i, j, k-1))*dzf_inv(k)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        return
    end subroutine calcVelGradTensor

    subroutine outVorticity(fn_prefix, nt)
        implicit none
        character(*), intent(in) :: fn_prefix
        integer, intent(in) :: nt

        real(fp) :: omega_x, omega_y, omega_z
        character(9) :: vartag = 'Vorticity'
        integer :: i, j, k

        if (.not. allocated(vor)) call allocVariable(nhalo_zero, sz, vartag, vor)

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k, omega_x, omega_y, omega_z)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            omega_x = vel_grad%wy(i, j, k) - vel_grad%vz(i, j, k)
            omega_y = vel_grad%uz(i, j, k) - vel_grad%wx(i, j, k)
            omega_z = vel_grad%vx(i, j, k) - vel_grad%uy(i, j, k)
            vor(i, j, k) = sqrt(omega_x**2 + omega_y**2 + omega_z**2)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        if (myrank == 0) write(*,'(A)') "PowerLLEL_Postprocess.NOTE.outVorticity: Writing inst. field <"//vartag//"> ..."
        call outputField(fn_prefix, nt, (/nx, ny, nz/), st, sz, nhalo_zero, vartag, vor)

        if (allocated(vor)) deallocate(vor)

        return
    end subroutine outVorticity

    subroutine outQ(fn_prefix, nt)
        implicit none
        character(*), intent(in) :: fn_prefix
        integer, intent(in) :: nt

        character(1) :: vartag = 'Q'
        integer :: i, j, k

        if (.not. allocated(q)) call allocVariable(nhalo_zero, sz, vartag, q)

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            q(i, j, k) = -(vel_grad%ux(i, j, k)**2 + vel_grad%vy(i, j, k)**2 + vel_grad%wz(i, j, k)**2) * 0.5 &
                         - vel_grad%uy(i, j, k) * vel_grad%vx(i, j, k) &
                         - vel_grad%uz(i, j, k) * vel_grad%wx(i, j, k) &
                         - vel_grad%vz(i, j, k) * vel_grad%wy(i, j, k)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        if (myrank == 0) write(*,'(A)') "PowerLLEL_Postprocess.NOTE.outQ: Writing inst. field <"//vartag//"> ..."
        call outputField(fn_prefix, nt, (/nx, ny, nz/), st, sz, nhalo_zero, vartag, q)

        if (allocated(q)) deallocate(q)

        return
    end subroutine outQ

    subroutine outLambda2(fn_prefix, nt)
        implicit none
        character(*), intent(in) :: fn_prefix
        integer, intent(in) :: nt

        character(7) :: vartag = 'Lambda2'
        real(fp) :: s11, s12, s13, s22, s23, s33
        real(fp) :: o12, o13, o23
        real(fp) :: a(3,3), x(3)
        integer :: i, j, k

        if (.not. allocated(lambda2)) call allocVariable(nhalo_zero, sz, vartag, lambda2)

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k, s11, s12, s13, s22, s23, s33, o12, o13, o23, a, x)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            ! Strain rate tensor S
            s11 = vel_grad%ux(i, j, k)
            s12 =(vel_grad%uy(i, j, k) + vel_grad%vx(i, j , k)) * 0.5
            s13 =(vel_grad%uz(i, j, k) + vel_grad%wx(i, j , k)) * 0.5
            s22 = vel_grad%vy(i, j, k)
            s23 =(vel_grad%vz(i, j, k) + vel_grad%wy(i, j , k)) * 0.5
            s33 = vel_grad%wz(i, j, k)
            ! Vorticity tensor \Omega
            o12 =(vel_grad%uy(i, j, k) - vel_grad%vx(i, j , k)) * 0.5
            o13 =(vel_grad%uz(i, j, k) - vel_grad%wx(i, j , k)) * 0.5
            o23 =(vel_grad%vz(i, j, k) - vel_grad%wy(i, j , k)) * 0.5

            ! S^2 + \Omega^2
            a(1,1) = s11*s11 + s12*s12 + s13*s13 - o12*o12 - o13*o13
            a(1,2) = s11*s12 + s12*s22 + s13*s23 - o13*o23
            a(1,3) = s11*s13 + s12*s23 + s13*s33 + o12*o23
            a(2,1) = a(1,2)
            a(2,2) = s12*s12 + s22*s22 + s23*s23 - o12*o12 - o23*o23
            a(2,3) = s12*s13 + s22*s23 + s23*s33 - o12*o13
            a(3,1) = a(1,3)
            a(3,2) = a(2,3)
            a(3,3) = s13*s13 + s23*s23 + s33*s33 - o13*o13 - o23*o23

            ! eigenvalue
            call eigenvalue(a, x)
            lambda2(i, j, k) = x(2)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        if (myrank == 0) write(*,'(A)') "PowerLLEL_Postprocess.NOTE.outLambda2: Writing inst. field <"//vartag//"> ..."
        call outputField(fn_prefix, nt, (/nx, ny, nz/), st, sz, nhalo_zero, vartag, lambda2)

        if (allocated(lambda2)) deallocate(lambda2)

        return
    end subroutine outLambda2

    ! gives back the real eigenvalues of a 3x3 real matrix
    subroutine eigenvalue(a, x)
        implicit none
        real(fp), intent(in ) :: a(3,3)
        real(fp), intent(out) :: x(3)

        real(fp) :: ca, cb, cc
        real(fp) :: Q, R
        complex(fp) :: AA, BB, cu

        cu = cmplx(0.0, 1.0)

        ca = -(a(1,1) + a(2,2) + a(3,3))
        cb = -(a(1,2)*a(2,1) + a(1,3)*a(3,1) + a(2,3)*a(3,2) &
             - a(1,1)*a(2,2) - a(1,1)*a(3,3) - a(2,2)*a(3,3))
        cc = -(a(1,2)*a(2,3)*a(3,1) - a(1,3)*a(2,2)*a(3,1) &
             + a(1,3)*a(2,1)*a(3,2) - a(1,1)*a(2,3)*a(3,2) &
             - a(1,2)*a(2,1)*a(3,3) + a(1,1)*a(2,2)*a(3,3))

        ! the solution to x^3 + A*x^2 + B*x + C = 0.

        Q = (ca*ca - 3.0*cb)/9.0
        R = (2.0*ca*ca*ca - 9.0*ca*cb + 27.0*cc)/54.0

        AA = -sign(1.0, R) * ( abs(R) + sqrt(cmplx(R*R-Q*Q*Q)) )**(1.0/3.0)

        if (AA == 0.0) then
            BB = 0.0
        else
            BB = cmplx(Q)/AA
        end if

        x(1) = (AA+BB) - ca/3.0
        x(2) = -0.5*(AA+BB) - ca/3.0 + cu*sqrt(3.0)/2.0*(AA-BB)
        x(3) = -0.5*(AA+BB) - ca/3.0 - cu*sqrt(3.0)/2.0*(AA-BB)

    end subroutine eigenvalue

end module mod_vortex_postproc