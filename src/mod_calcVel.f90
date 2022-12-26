#define MEAN Mean_h

module mod_calcVel
    use mod_type,       only: fp
    use mod_parameters, only: re_inv, dt, is_forced, vel_force, nx, ny, nz, lz, nhalo, u_crf, &
                              smooth_wall_visc
    use mod_mpi,        only: sz
    use mod_mesh,       only: dx_inv, dy_inv, dzc_inv, dzf, dzf_inv, dzflzi, dzf_global, dzc, visc_dzf_inv
    use mod_utils,      only: MEAN
#ifdef NB_HALO
    use mod_updateBound,only: imposeBCVel, updateHaloISend, updateHaloIRecv, updateHaloWaitall
#endif
    !$ use omp_lib
#ifdef GPTL
    use gptl
#endif

    implicit none

#ifdef NB_HALO
    include 'mpif.h'
#endif

contains

    subroutine transform2CRF(vel_crf, nhalo, sz, vel, vel_force)
        implicit none
        real(fp), intent(in) :: vel_crf
        integer, dimension(6), intent(in) :: nhalo
        integer, dimension(3), intent(in) :: sz
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: vel
        real(fp), intent(inout) :: vel_force

        integer :: i, j, k

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            vel(i, j, k) = vel(i, j, k) - vel_crf
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        vel_force = vel_force - vel_crf

        return
    end subroutine transform2CRF

    subroutine timeIntVelRK1_kernel(st, en, u, v, w, u1, v1, w1, u_crf)
        implicit none
        integer,  dimension(3), intent(in) :: st, en
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in), contiguous :: u, v, w
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(out), contiguous :: u1, v1, w1
        real(fp), intent(in) :: u_crf

        real(fp) :: conv, visc
        real(fp) :: uw, ue, us, un, ub, ut
        real(fp) :: ww, we, ws, wn, wb, wt
        real(fp) :: vw, ve, vs, vn, vb, vt
        real(fp) :: dqw, dqe, dqs, dqn, dqb, dqt
        integer  :: i, j, k

        !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2) &
        !$OMP DEFAULT(SHARED) &
        !$OMP PRIVATE(i,j,k,ue,uw,un,us,ut,ub,ve,vw,vn,vs,vt,vb,we,ww,ws,wn,wt,wb, &
        !$OMP         dqe,dqw,dqn,dqs,dqt,dqb,conv,visc)
        do k = st(3), en(3)
        do j = st(2), en(2)
        do i = st(1), en(1)

            ue  = (u(i  ,j  ,k  )+u(i+1,j  ,k  ))
            uw  = (u(i-1,j  ,k  )+u(i  ,j  ,k  ))
            un  = (u(i  ,j  ,k  )+u(i  ,j+1,k  ))
            us  = (u(i  ,j-1,k  )+u(i  ,j  ,k  ))
            ut  = (u(i  ,j  ,k  )*dzf(k+1) + u(i,j,k+1)*dzf(k))*dzc_inv(k  )
            ub  = (u(i  ,j  ,k  )*dzf(k-1) + u(i,j,k-1)*dzf(k))*dzc_inv(k-1)
            vn  = (v(i  ,j  ,k  )+v(i+1,j  ,k  ))
            vs  = (v(i  ,j-1,k  )+v(i+1,j-1,k  ))
            wt  = (w(i  ,j  ,k  )+w(i+1,j  ,k  ))
            wb  = (w(i  ,j  ,k-1)+w(i+1,j  ,k-1))
            dqe = (u(i+1,j  ,k  )-u(i  ,j  ,k  ))*dx_inv
            dqw = (u(i  ,j  ,k  )-u(i-1,j  ,k  ))*dx_inv
            dqn = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))*dy_inv
            dqs = (u(i  ,j  ,k  )-u(i  ,j-1,k  ))*dy_inv
            dqt = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))*dzc_inv(k  )
            dqb = (u(i  ,j  ,k  )-u(i  ,j  ,k-1))*dzc_inv(k-1)
            conv = 0.25_fp*( (ue*ue-uw*uw)*dx_inv + (un*vn-us*vs)*dy_inv + (ut*wt-ub*wb)*dzf_inv(k) )
            ! add a term induced by the convecting reference frame
            conv = conv + 0.5_fp*u_crf*(ue-uw)*dx_inv
            visc = ((dqe-dqw)*dx_inv + (dqn-dqs)*dy_inv + (dqt-dqb)*dzf_inv(k))*re_inv
            u1(i,j,k) = u(i,j,k) + dt*(visc - conv)

            ve  = (v(i  ,j  ,k  )+v(i+1,j  ,k  ))
            vw  = (v(i-1,j  ,k  )+v(i  ,j  ,k  ))
            vn  = (v(i  ,j  ,k  )+v(i  ,j+1,k  ))
            vs  = (v(i  ,j-1,k  )+v(i  ,j  ,k  ))
            vt  = (v(i  ,j  ,k  )*dzf(k+1) + v(i,j,k+1)*dzf(k))*dzc_inv(k  )
            vb  = (v(i  ,j  ,k  )*dzf(k-1) + v(i,j,k-1)*dzf(k))*dzc_inv(k-1)
            ue  = (u(i  ,j  ,k  )+u(i  ,j+1,k  ))
            uw  = (u(i-1,j  ,k  )+u(i-1,j+1,k  ))
            wt  = (w(i  ,j  ,k  )+w(i  ,j+1,k  ))
            wb  = (w(i  ,j  ,k-1)+w(i  ,j+1,k-1))
            dqe = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))*dx_inv
            dqw = (v(i  ,j  ,k  )-v(i-1,j  ,k  ))*dx_inv
            dqn = (v(i  ,j+1,k  )-v(i  ,j  ,k  ))*dy_inv
            dqs = (v(i  ,j  ,k  )-v(i  ,j-1,k  ))*dy_inv
            dqt = (v(i  ,j  ,k+1)-v(i  ,j  ,k  ))*dzc_inv(k  )
            dqb = (v(i  ,j  ,k  )-v(i  ,j  ,k-1))*dzc_inv(k-1)
            conv = 0.25_fp*( (ue*ve-uw*vw)*dx_inv + (vn*vn-vs*vs)*dy_inv + (wt*vt-wb*vb)*dzf_inv(k) )
            ! add a term induced by the convecting reference frame
            conv = conv + 0.5_fp*u_crf*(ve-vw)*dx_inv
            visc = ((dqe-dqw)*dx_inv + (dqn-dqs)*dy_inv + (dqt-dqb)*dzf_inv(k))*re_inv
            v1(i,j,k) = v(i,j,k) + dt*(visc - conv)

            we  = (w(i  ,j  ,k  )+w(i+1,j  ,k  ))
            ww  = (w(i-1,j  ,k  )+w(i  ,j  ,k  ))
            wn  = (w(i  ,j  ,k  )+w(i  ,j+1,k  ))
            ws  = (w(i  ,j-1,k  )+w(i  ,j  ,k  ))
            wt  = (w(i  ,j  ,k  )+w(i  ,j  ,k+1))
            wb  = (w(i  ,j  ,k  )+w(i  ,j  ,k-1))
            ue  = (u(i  ,j  ,k  )*dzf(k+1) + u(i  ,j  ,k+1)*dzf(k))*dzc_inv(k  )
            uw  = (u(i-1,j  ,k  )*dzf(k+1) + u(i-1,j  ,k+1)*dzf(k))*dzc_inv(k  )
            vn  = (v(i  ,j  ,k  )*dzf(k+1) + v(i  ,j  ,k+1)*dzf(k))*dzc_inv(k  )
            vs  = (v(i  ,j-1,k  )*dzf(k+1) + v(i  ,j-1,k+1)*dzf(k))*dzc_inv(k  )
            dqe = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))*dx_inv
            dqw = (w(i  ,j  ,k  )-w(i-1,j  ,k  ))*dx_inv
            dqn = (w(i  ,j+1,k  )-w(i  ,j  ,k  ))*dy_inv
            dqs = (w(i  ,j  ,k  )-w(i  ,j-1,k  ))*dy_inv
            dqt = (w(i  ,j  ,k+1)-w(i  ,j  ,k  ))*dzf_inv(k+1)
            dqb = (w(i  ,j  ,k  )-w(i  ,j  ,k-1))*dzf_inv(k)
            conv = 0.25_fp*( (ue*we-uw*ww)*dx_inv + (vn*wn-vs*ws)*dy_inv + (wt*wt-wb*wb)*dzc_inv(k) )
            ! add a term induced by the convecting reference frame
            conv = conv + 0.5_fp*u_crf*(we-ww)*dx_inv
            visc = ((dqe-dqw)*dx_inv + (dqn-dqs)*dy_inv + (dqt-dqb)*dzc_inv(k))*re_inv
            w1(i,j,k) = w(i,j,k) + dt*(visc - conv)

        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        return
    end subroutine timeIntVelRK1_kernel

    subroutine timeIntVelRK1(u, v, w, u1, v1, w1, u_crf)
        implicit none
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in), contiguous :: u, v, w
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(out), contiguous :: u1, v1, w1
        real(fp), intent(in) :: u_crf
        
#ifndef NB_HALO
        call timeIntVelRK1_kernel((/1,1,1/), sz, u, v, w, u1, v1, w1, u_crf)
#else
        integer :: tag_u = 1, tag_v = 2, tag_w = 3
        integer, dimension(8) :: isend_req_u, irecv_req_u
        integer, dimension(8) :: isend_req_v, irecv_req_v
        integer, dimension(8) :: isend_req_w, irecv_req_w

        integer :: ist, jst, kst
        integer :: y0_send_st, y0_send_en
        integer :: y1_send_st, y1_send_en
        integer :: z0_send_st, z0_send_en
        integer :: z1_send_st, z1_send_en

        integer :: ret

        ist = 1-nhalo(1)
        jst = 1-nhalo(3)
        kst = 1-nhalo(5)

        z0_send_st = 1;         z0_send_en = nhalo(6)
        z1_send_st = sz(3)+kst; z1_send_en = sz(3)
        y0_send_st = 1;         y0_send_en = nhalo(4)
        y1_send_st = sz(2)+jst; y1_send_en = sz(2)

#ifdef GPTL
        ret = gptlstart('--Update halo vel')
#endif

        call updateHaloIRecv(nhalo, tag_u, u1, irecv_req_u)
        call updateHaloIRecv(nhalo, tag_v, v1, irecv_req_v)
        call updateHaloIRecv(nhalo, tag_w, w1, irecv_req_w)

#ifdef GPTL
        ret = gptlstop('--Update halo vel')
        ret = gptlstart('--uvw1 comp')
#endif
        ! *** bottom/top ***
        call timeIntVelRK1_kernel((/1,1,z0_send_st/), (/sz(1),sz(2),z0_send_en/), &
                                  u, v, w, u1, v1, w1, u_crf)
        call timeIntVelRK1_kernel((/1,1,max(z1_send_st,z0_send_en+1)/), (/sz(1),sz(2),z1_send_en/), &
                                  u, v, w, u1, v1, w1, u_crf)
        ! *** south/north ***
        call timeIntVelRK1_kernel((/1,y0_send_st,z0_send_en+1/), (/sz(1),y0_send_en,z1_send_st-1/), &
                                  u, v, w, u1, v1, w1, u_crf)
        call timeIntVelRK1_kernel((/1,max(y1_send_st,y0_send_en+1),z0_send_en+1/), (/sz(1),y1_send_en,z1_send_st-1/), &
                                  u, v, w, u1, v1, w1, u_crf)
#ifdef GPTL
        ret = gptlstop('--uvw1 comp')
        ret = gptlstart('--Update halo vel')
#endif

        call updateHaloISend(nhalo, tag_u, u1, isend_req_u)
        call updateHaloISend(nhalo, tag_v, v1, isend_req_v)
        call updateHaloISend(nhalo, tag_w, w1, isend_req_w)

#ifdef GPTL
        ret = gptlstop('--Update halo vel')
        ret = gptlstart('--uvw1 comp')
#endif
        ! *** inner region ***
        call timeIntVelRK1_kernel((/1,nhalo(4)+1,nhalo(6)+1/), (/sz(1),sz(2)-nhalo(3),sz(3)-nhalo(5)/), &
                                  u, v, w, u1, v1, w1, u_crf)
#ifdef GPTL
        ret = gptlstop('--uvw1 comp')
        ret = gptlstart('--Update halo vel')
#endif

        call updateHaloWaitall(isend_req_u, irecv_req_u)
        call updateHaloWaitall(isend_req_v, irecv_req_v)
        call updateHaloWaitall(isend_req_w, irecv_req_w)

#ifdef GPTL
        ret = gptlstop('--Update halo vel')
        ret = gptlstart('--Impose BC vel')
#endif
        call imposeBCVel(u1, v1, w1, u_crf)
#ifdef GPTL
        ret = gptlstop('--Impose BC vel')
#endif

#endif

        return
    end subroutine timeIntVelRK1
    
    subroutine timeIntVelRK2_kernel(st, en, u, v, w, u1, v1, w1, u_crf)
        implicit none
        integer,  dimension(3), intent(in) :: st, en
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout), contiguous :: u, v, w
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in), contiguous :: u1, v1, w1
        real(fp), intent(in) :: u_crf

        real(fp) :: conv, visc
        real(fp) :: uw, ue, us, un, ub, ut
        real(fp) :: ww, we, ws, wn, wb, wt
        real(fp) :: vw, ve, vs, vn, vb, vt
        real(fp) :: dqw, dqe, dqs, dqn, dqb, dqt
        integer  :: i, j, k

        !$OMP PARALLEL DO SCHEDULE(STATIC) COLLAPSE(2) &
        !$OMP DEFAULT(SHARED) &
        !$OMP PRIVATE(i,j,k,ue,uw,un,us,ut,ub,ve,vw,vn,vs,vt,vb,we,ww,ws,wn,wt,wb, &
        !$OMP         dqe,dqw,dqn,dqs,dqt,dqb,conv,visc)
        do k = st(3), en(3)
        do j = st(2), en(2)
        do i = st(1), en(1)
        
            ue  = (u1(i  ,j  ,k  )+u1(i+1,j  ,k  ))
            uw  = (u1(i-1,j  ,k  )+u1(i  ,j  ,k  ))
            un  = (u1(i  ,j  ,k  )+u1(i  ,j+1,k  ))
            us  = (u1(i  ,j-1,k  )+u1(i  ,j  ,k  ))
            ut  = (u1(i  ,j  ,k  )*dzf(k+1) + u1(i,j,k+1)*dzf(k))*dzc_inv(k  )
            ub  = (u1(i  ,j  ,k  )*dzf(k-1) + u1(i,j,k-1)*dzf(k))*dzc_inv(k-1)
            vn  = (v1(i  ,j  ,k  )+v1(i+1,j  ,k  ))
            vs  = (v1(i  ,j-1,k  )+v1(i+1,j-1,k  ))
            wt  = (w1(i  ,j  ,k  )+w1(i+1,j  ,k  ))
            wb  = (w1(i  ,j  ,k-1)+w1(i+1,j  ,k-1))
            dqe = (u1(i+1,j  ,k  )-u1(i  ,j  ,k  ))*dx_inv
            dqw = (u1(i  ,j  ,k  )-u1(i-1,j  ,k  ))*dx_inv
            dqn = (u1(i  ,j+1,k  )-u1(i  ,j  ,k  ))*dy_inv
            dqs = (u1(i  ,j  ,k  )-u1(i  ,j-1,k  ))*dy_inv
            dqt = (u1(i  ,j  ,k+1)-u1(i  ,j  ,k  ))*dzc_inv(k  )
            dqb = (u1(i  ,j  ,k  )-u1(i  ,j  ,k-1))*dzc_inv(k-1)
            conv = 0.25_fp*( (ue*ue-uw*uw)*dx_inv + (un*vn-us*vs)*dy_inv + (ut*wt-ub*wb)*dzf_inv(k) )
            ! add a term induced by the convecting reference frame
            conv = conv + 0.5_fp*u_crf*(ue-uw)*dx_inv
            visc = ((dqe-dqw)*dx_inv + (dqn-dqs)*dy_inv + (dqt-dqb)*dzf_inv(k))*re_inv
            u(i,j,k) = (u1(i,j,k) + dt*(visc-conv) + u(i,j,k))*0.5_fp

            ve  = (v1(i  ,j  ,k  )+v1(i+1,j  ,k  ))
            vw  = (v1(i-1,j  ,k  )+v1(i  ,j  ,k  ))
            vn  = (v1(i  ,j  ,k  )+v1(i  ,j+1,k  ))
            vs  = (v1(i  ,j-1,k  )+v1(i  ,j  ,k  ))
            vt  = (v1(i  ,j  ,k  )*dzf(k+1) + v1(i,j,k+1)*dzf(k))*dzc_inv(k  )
            vb  = (v1(i  ,j  ,k  )*dzf(k-1) + v1(i,j,k-1)*dzf(k))*dzc_inv(k-1)
            ue  = (u1(i  ,j  ,k  )+u1(i  ,j+1,k  ))
            uw  = (u1(i-1,j  ,k  )+u1(i-1,j+1,k  ))
            wt  = (w1(i  ,j  ,k  )+w1(i  ,j+1,k  ))
            wb  = (w1(i  ,j  ,k-1)+w1(i  ,j+1,k-1))
            dqe = (v1(i+1,j  ,k  )-v1(i  ,j  ,k  ))*dx_inv
            dqw = (v1(i  ,j  ,k  )-v1(i-1,j  ,k  ))*dx_inv
            dqn = (v1(i  ,j+1,k  )-v1(i  ,j  ,k  ))*dy_inv
            dqs = (v1(i  ,j  ,k  )-v1(i  ,j-1,k  ))*dy_inv
            dqt = (v1(i  ,j  ,k+1)-v1(i  ,j  ,k  ))*dzc_inv(k  )
            dqb = (v1(i  ,j  ,k  )-v1(i  ,j  ,k-1))*dzc_inv(k-1)
            conv = 0.25_fp*( (ue*ve-uw*vw)*dx_inv + (vn*vn-vs*vs)*dy_inv + (wt*vt-wb*vb)*dzf_inv(k) )
            ! add a term induced by the convecting reference frame
            conv = conv + 0.5_fp*u_crf*(ve-vw)*dx_inv
            visc = ((dqe-dqw)*dx_inv + (dqn-dqs)*dy_inv + (dqt-dqb)*dzf_inv(k))*re_inv
            v(i,j,k) = (v1(i,j,k) + dt*(visc-conv) + v(i,j,k))*0.5_fp

            we  = (w1(i  ,j  ,k  )+w1(i+1,j  ,k  ))
            ww  = (w1(i-1,j  ,k  )+w1(i  ,j  ,k  ))
            wn  = (w1(i  ,j  ,k  )+w1(i  ,j+1,k  ))
            ws  = (w1(i  ,j-1,k  )+w1(i  ,j  ,k  ))
            wt  = (w1(i  ,j  ,k  )+w1(i  ,j  ,k+1))
            wb  = (w1(i  ,j  ,k  )+w1(i  ,j  ,k-1))
            ue  = (u1(i  ,j  ,k  )*dzf(k+1) + u1(i  ,j  ,k+1)*dzf(k))*dzc_inv(k  )
            uw  = (u1(i-1,j  ,k  )*dzf(k+1) + u1(i-1,j  ,k+1)*dzf(k))*dzc_inv(k  )
            vn  = (v1(i  ,j  ,k  )*dzf(k+1) + v1(i  ,j  ,k+1)*dzf(k))*dzc_inv(k  )
            vs  = (v1(i  ,j-1,k  )*dzf(k+1) + v1(i  ,j-1,k+1)*dzf(k))*dzc_inv(k  )
            dqe = (w1(i+1,j  ,k  )-w1(i  ,j  ,k  ))*dx_inv
            dqw = (w1(i  ,j  ,k  )-w1(i-1,j  ,k  ))*dx_inv
            dqn = (w1(i  ,j+1,k  )-w1(i  ,j  ,k  ))*dy_inv
            dqs = (w1(i  ,j  ,k  )-w1(i  ,j-1,k  ))*dy_inv
            dqt = (w1(i  ,j  ,k+1)-w1(i  ,j  ,k  ))*dzf_inv(k+1)
            dqb = (w1(i  ,j  ,k  )-w1(i  ,j  ,k-1))*dzf_inv(k)
            conv = 0.25_fp*( (ue*we-uw*ww)*dx_inv + (vn*wn-vs*ws)*dy_inv + (wt*wt-wb*wb)*dzc_inv(k) )
            ! add a term induced by the convecting reference frame
            conv = conv + 0.5_fp*u_crf*(we-ww)*dx_inv
            visc = ((dqe-dqw)*dx_inv + (dqn-dqs)*dy_inv + (dqt-dqb)*dzc_inv(k))*re_inv
            w(i,j,k) = (w1(i,j,k) + dt*(visc-conv) + w(i,j,k))*0.5_fp

        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        return
    end subroutine timeIntVelRK2_kernel

    subroutine timeIntVelRK2(u, v, w, u1, v1, w1, u_crf)
        implicit none
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout), contiguous :: u, v, w
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(in), contiguous :: u1, v1, w1
        real(fp), intent(in) :: u_crf

#ifndef NB_HALO
        call timeIntVelRK2_kernel((/1,1,1/), sz, u, v, w, u1, v1, w1, u_crf)
#else
        integer :: tag_u = 1, tag_v = 2, tag_w = 3
        integer, dimension(8) :: isend_req_u, irecv_req_u
        integer, dimension(8) :: isend_req_v, irecv_req_v
        integer, dimension(8) :: isend_req_w, irecv_req_w

        integer :: ist, jst, kst
        integer :: y0_send_st, y0_send_en
        integer :: y1_send_st, y1_send_en
        integer :: z0_send_st, z0_send_en
        integer :: z1_send_st, z1_send_en

        integer :: ret

        ist = 1-nhalo(1)
        jst = 1-nhalo(3)
        kst = 1-nhalo(5)

        z0_send_st = 1;         z0_send_en = nhalo(6)
        z1_send_st = sz(3)+kst; z1_send_en = sz(3)
        y0_send_st = 1;         y0_send_en = nhalo(4)
        y1_send_st = sz(2)+jst; y1_send_en = sz(2)

#ifdef GPTL
        ret = gptlstart('--Update halo vel')
#endif

        call updateHaloIRecv(nhalo, tag_u, u, irecv_req_u)
        call updateHaloIRecv(nhalo, tag_v, v, irecv_req_v)
        call updateHaloIRecv(nhalo, tag_w, w, irecv_req_w)

#ifdef GPTL
        ret = gptlstop('--Update halo vel')
        ret = gptlstart('--uvw2 comp')
#endif
        ! *** bottom/top ***
        call timeIntVelRK2_kernel((/1,1,z0_send_st/), (/sz(1),sz(2),z0_send_en/), &
                                  u, v, w, u1, v1, w1, u_crf)
        call timeIntVelRK2_kernel((/1,1,max(z1_send_st,z0_send_en+1)/), (/sz(1),sz(2),z1_send_en/), &
                                  u, v, w, u1, v1, w1, u_crf)
        ! *** south/north ***
        call timeIntVelRK2_kernel((/1,y0_send_st,z0_send_en+1/), (/sz(1),y0_send_en,z1_send_st-1/), &
                                  u, v, w, u1, v1, w1, u_crf)
        call timeIntVelRK2_kernel((/1,max(y1_send_st,y0_send_en+1),z0_send_en+1/), (/sz(1),y1_send_en,z1_send_st-1/), &
                                  u, v, w, u1, v1, w1, u_crf)
#ifdef GPTL
        ret = gptlstop('--uvw2 comp')
        ret = gptlstart('--Update halo vel')
#endif

        call updateHaloISend(nhalo, tag_u, u, isend_req_u)
        call updateHaloISend(nhalo, tag_v, v, isend_req_v)
        call updateHaloISend(nhalo, tag_w, w, isend_req_w)

#ifdef GPTL
        ret = gptlstop('--Update halo vel')
        ret = gptlstart('--uvw2 comp')
#endif
        ! *** inner region ***
        call timeIntVelRK2_kernel((/1,nhalo(4)+1,nhalo(6)+1/), (/sz(1),sz(2)-nhalo(3),sz(3)-nhalo(5)/), &
                                  u, v, w, u1, v1, w1, u_crf)
#ifdef GPTL
        ret = gptlstop('--uvw2 comp')
        ret = gptlstart('--Update halo vel')
#endif

        call updateHaloWaitall(isend_req_u, irecv_req_u)
        call updateHaloWaitall(isend_req_v, irecv_req_v)
        call updateHaloWaitall(isend_req_w, irecv_req_w)

#ifdef GPTL
        ret = gptlstop('--Update halo vel')
        ret = gptlstart('--Impose BC vel')
#endif
        call imposeBCVel(u, v, w, u_crf)
#ifdef GPTL
        ret = gptlstop('--Impose BC vel')
#endif

#endif

        return
    end subroutine timeIntVelRK2

    subroutine correctVel(p, u, v, w)
        implicit none
        real(fp), dimension(0:,0:,0:), intent(in   ) :: p
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: u, v, w

        real(fp) :: dtdxi, dtdyi
        integer :: i, j, k

        dtdxi = dt*dx_inv
        dtdyi = dt*dy_inv

        !$OMP PARALLEL DO SCHEDULE(STATIC) &
        !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
        do k = 1, sz(3)
        do j = 1, sz(2)
        do i = 1, sz(1)
            u(i, j, k) = u(i, j, k) - (p(i+1, j, k)-p(i, j, k))*dtdxi
            v(i, j, k) = v(i, j, k) - (p(i, j+1, k)-p(i, j, k))*dtdyi
            w(i, j, k) = w(i, j, k) - (p(i, j, k+1)-p(i, j, k))*dt*dzc_inv(k)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO

        return
    end subroutine correctVel

    subroutine forceVel(u, v, w)
        implicit none
        real(fp), dimension(1-nhalo(1):,1-nhalo(3):,1-nhalo(5):), intent(inout) :: u, v, w

        real(fp) :: vel_mean_x, vel_mean_y, vel_mean_z
        real(fp) :: force_x, force_y, force_z
        integer :: i, j, k

        if (is_forced(1)) then
            vel_mean_x = MEAN(nhalo, sz, nx, ny, dzflzi, u)
            force_x = vel_force(1) - vel_mean_x
            !$OMP PARALLEL DO SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
            do k = 1, sz(3)
            do j = 1, sz(2)
            do i = 1, sz(1)
                u(i, j, k) = u(i, j, k) + force_x
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif

        if (is_forced(2)) then
            vel_mean_y = MEAN(nhalo, sz, nx, ny, dzflzi, v)
            force_y = vel_force(2) - vel_mean_y
            !$OMP PARALLEL DO SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
            do k = 1, sz(3)
            do j = 1, sz(2)
            do i = 1, sz(1)
                v(i, j, k) = v(i, j, k) + force_y
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif

        if (is_forced(3)) then
            vel_mean_z = MEAN(nhalo, sz, nx, ny, dzflzi, w)
            force_z = vel_force(3) - vel_mean_z
            !$OMP PARALLEL DO SCHEDULE(STATIC) &
            !$OMP DEFAULT(SHARED) PRIVATE(i, j, k)
            do k = 1, sz(3)
            do j = 1, sz(2)
            do i = 1, sz(1)
                w(i, j, k) = w(i, j, k) + force_z
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif

        return
    end subroutine forceVel

end module mod_calcVel