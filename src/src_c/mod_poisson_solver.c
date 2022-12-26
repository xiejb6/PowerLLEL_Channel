#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#ifdef GPTL
#include <gptl.h>
#endif
#include "decomp_2d.h"
#include "mod_fft.h"
#include "memory.h"

#ifdef USE_MKL
    fft_mkl_plan_ptr plan[2][2];
#else
    fftw_plan plan[2][2];
#endif
double fft_normfactor;
int *xsize, *ysize, *zsize;
double *var_xpen, *var_ypen, *var_zpen;
double *a, *b, *c;
int *sz_trid, *st_trid, *neighbor_trid;

#ifdef _PDD
MPI_Comm COMM_CART_PDD;
double *w_pdd, *v_pdd, *tmp_v_pdd;
double *y1_pdd, *y2_pdd, *y3_pdd;
double *tmp_var_pdd;
size_t tmp_var_pdd_len;
#endif

size_t var_offset1;
size_t var_offset2;
size_t var_xpen_offset1;
size_t var_xpen_offset2;
size_t var_ypen_offset1;
size_t var_ypen_offset2;

static void set_trid_coeff(int n, double *dzf, char bctype[], char c_or_f, int neighbor[2], 
                           double a[n], double b[n], double c[n]) {
    int dzf_st = -1;
    if (c_or_f == 'c') {
        for (int k = 0; k < n; k++) {
            a[k] = 2.0 / (dzf[k - dzf_st] * (dzf[k-1 - dzf_st] + dzf[k - dzf_st]));
            c[k] = 2.0 / (dzf[k - dzf_st] * (dzf[k+1 - dzf_st] + dzf[k - dzf_st]));
        }
    } else if (c_or_f == 'f') {
        for (int k = 0; k < n; k++) {
            a[k] = 2.0 / (dzf[k   - dzf_st] * (dzf[k+1 - dzf_st] + dzf[k - dzf_st]));
            c[k] = 2.0 / (dzf[k+1 - dzf_st] * (dzf[k+1 - dzf_st] + dzf[k - dzf_st]));
        }
    }

    for (int i = 0; i < n; i++) {
        b[i] = - a[i] - c[i];
    }

    // coefficients correction according to BC types
    double factor;
    char bc;
    if (neighbor[0] == MPI_PROC_NULL) {
        bc = bctype[0];
        if (bc == 'P') {
            factor = 0.0;
        } else if (bc == 'D') {
            factor = -1.0;
        } else if (bc == 'N') {
            factor = 1.0;
        }

        if (c_or_f == 'c') {
            b[0] += factor * a[0];
            a[0] *= fabs(factor) - 1.0;
        } else if (c_or_f == 'f') {
            if (bc == 'N') {
                b[0] += factor * a[0];
                a[0] *= fabs(factor) - 1.0;
            }
        }
    }
    
    if (neighbor[1] == MPI_PROC_NULL) {
        bc = bctype[1];
        if (bc == 'P') {
            factor = 0.0;
        } else if (bc == 'D') {
            factor = -1.0;
        } else if (bc == 'N') {
            factor = 1.0;
        }

        if (c_or_f == 'c') {
            b[n-1] += factor * c[n-1];
            c[n-1] *= fabs(factor) - 1.0;
        } else if (c_or_f == 'f') {
            if (bc == 'N') {
                b[n-1] += factor * c[n-1];
                c[n-1] *= fabs(factor) - 1.0;
            }
        }
    }
}

#ifdef _PDD
static void init_pdd_array(int sz[3]) {
    v_pdd = (double *)aligned_malloc(sizeof(double) * sz[0] * sz[1] * sz[2], MEM_ALIGN_SIZE);
    w_pdd = (double *)aligned_malloc(sizeof(double) * sz[0] * sz[1] * sz[2], MEM_ALIGN_SIZE);
    tmp_var_pdd_len = sz[0]*sz[1];
    y1_pdd = (double *)aligned_malloc(sizeof(double) * tmp_var_pdd_len, MEM_ALIGN_SIZE);
    y2_pdd = (double *)aligned_malloc(sizeof(double) * tmp_var_pdd_len, MEM_ALIGN_SIZE);
    y3_pdd = (double *)aligned_malloc(sizeof(double) * tmp_var_pdd_len, MEM_ALIGN_SIZE);
    tmp_v_pdd = (double *)aligned_malloc(sizeof(double) * tmp_var_pdd_len, MEM_ALIGN_SIZE);
    tmp_var_pdd = (double *)aligned_malloc(sizeof(double) * tmp_var_pdd_len, MEM_ALIGN_SIZE);

    memset(v_pdd, 0, sizeof(double) * sz[0] * sz[1] * sz[2]);
    memset(w_pdd, 0, sizeof(double) * sz[0] * sz[1] * sz[2]);
    memset(y1_pdd, 0, sizeof(double) * tmp_var_pdd_len);
    memset(y2_pdd, 0, sizeof(double) * tmp_var_pdd_len);
    memset(y3_pdd, 0, sizeof(double) * tmp_var_pdd_len);

    int ie = sz[0];
    int je = sz[1];
    int ke = sz[2];
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < je; j++) {
        for (int i = 0; i < ie; i++) {
            v_pdd[i + j * ie +      0 * ie * je] = a[0];
            w_pdd[i + j * ie + (ke-1) * ie * je] = c[ke-1];
        }
    }
    a[0] = 0.0;
    c[ke-1] = 0.0;

    #pragma omp parallel
    {
        for (int k = 1; k < ke; k++) {
            #pragma omp for schedule(static)
            for (int j = 0; j < je; j++) {
                for (int i = 0; i < ie; i++) {
                    double a_tmp = a[k] / b[i + j * ie + (k-1) * ie * je];
                    v_pdd[i + j * ie + k * ie * je] -= a_tmp * v_pdd[i + j * ie + (k-1) * ie * je];
                    w_pdd[i + j * ie + k * ie * je] -= a_tmp * w_pdd[i + j * ie + (k-1) * ie * je];
                }
            }
        }

        #pragma omp for schedule(static)
        for (int j = 0; j < je; j++) {
            for (int i = 0; i < ie; i++) {
                if (b[i + j * ie + (ke-1) * ie * je] != 0.0) {
                    v_pdd[i + j * ie + (ke-1) * ie * je] /= b[i + j * ie + (ke-1) * ie * je];
                    w_pdd[i + j * ie + (ke-1) * ie * je] /= b[i + j * ie + (ke-1) * ie * je];
                } else {
                    v_pdd[i + j * ie + (ke-1) * ie * je] = 0.0;
                    w_pdd[i + j * ie + (ke-1) * ie * je] = 0.0;
                }
            }
        }

        for (int k = ke-2; k >= 0; k--) {
            #pragma omp for schedule(static)
            for (int j = 0; j < je; j++) {
                for (int i = 0; i < ie; i++) {
                    v_pdd[i + j * ie + k * ie * je] = (v_pdd[i + j * ie + k * ie * je] - c[k] * v_pdd[i + j * ie + (k+1) * ie * je]) / b[i + j * ie + k * ie * je];
                    w_pdd[i + j * ie + k * ie * je] = (w_pdd[i + j * ie + k * ie * je] - c[k] * w_pdd[i + j * ie + (k+1) * ie * je]) / b[i + j * ie + k * ie * je];
                }
            }
        }
    }

    MPI_Sendrecv(v_pdd, tmp_var_pdd_len, MPI_DOUBLE, neighbor_trid[4], 100,
                 tmp_v_pdd, tmp_var_pdd_len, MPI_DOUBLE, neighbor_trid[5], 100,
                 COMM_CART_PDD, MPI_STATUS_IGNORE);
}
#endif

void init_poisson_solver(int nx_global, int ny_global, int nz_global, 
                         double dx, double dy, double *dzf_global,
                         char bctype_x[], char bctype_y[], char bctype_z[], 
                         int neighbor_xyz[3][6]) {

    xsize = decomp_main.xsz;
    ysize = decomp_main.ysz;
    zsize = decomp_main.zsz;

    var_offset1 = (xsize[0]+2) * (xsize[1]+2);
    var_offset2 = (xsize[0]+2);
    var_xpen_offset1 = xsize[0] * xsize[1];
    var_xpen_offset2 = xsize[0];
    var_ypen_offset1 = ysize[0] * ysize[1];
    var_ypen_offset2 = ysize[0];

    var_xpen = (double*)aligned_malloc(sizeof(double)*xsize[0]*xsize[1]*xsize[2] + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
    if (xsize[0] == ysize[0] && xsize[1] == ysize[1] && xsize[2] == ysize[2]) {
        var_ypen = var_xpen;
    }
    else {
        var_ypen = (double*)aligned_malloc(sizeof(double)*ysize[0]*ysize[1]*ysize[2] + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
    }
    if (ysize[0] == zsize[0] && ysize[1] == zsize[1] && ysize[2] == zsize[2]) {
        var_zpen = var_ypen;
    }
    else {
        var_zpen = (double*)aligned_malloc(sizeof(double)*zsize[0]*zsize[1]*zsize[2] + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
    }

    // determine a decomposition mode used for solving the tridiagonal system
#ifdef _PDD
    sz_trid = decomp_main.ysz;
    st_trid = decomp_main.yst;
    neighbor_trid = &neighbor_xyz[1][0];
#else
    sz_trid = decomp_main.zsz;
    st_trid = decomp_main.zst;
    neighbor_trid = &neighbor_xyz[2][0];
#endif

    // initialize FFT
#ifdef USE_MKL
    init_fft(xsize, ysize, bctype_x, bctype_y, plan, &fft_normfactor);
#else
    init_fft(xsize, ysize, bctype_x, bctype_y, var_xpen, var_ypen, plan, &fft_normfactor);
#endif

    // calculate eigenvalues corresponding to BC types
    double *lambdax = (double *)malloc(sz_trid[0] * sizeof(double));
    double *lambday = (double *)malloc(sz_trid[1] * sizeof(double));
    double *lambdaxy = (double *)malloc(sz_trid[0] * sz_trid[1] * sizeof(double));
    get_eigen_values(st_trid[0], sz_trid[0], nx_global, bctype_x, lambdax);
    for (int i = 0; i < sz_trid[0]; i++) {
        lambdax[i] /= (dx*dx);
    }
    get_eigen_values(st_trid[1], sz_trid[1], ny_global, bctype_y, lambday);
    for (int i = 0; i < sz_trid[1]; i++) {
        lambday[i] /= (dy*dy);
    }
    for (int j = 0; j < sz_trid[1]; j++) {
        for (int i = 0; i < sz_trid[0]; i++) {
            lambdaxy[i + j * sz_trid[0]] = lambdax[i] + lambday[j];
        }
    }
    
    // calculate coefficients of tridiagonal systems
    a = (double *)aligned_malloc(sizeof(double)*sz_trid[2] + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
    c = (double *)aligned_malloc(sizeof(double)*sz_trid[2] + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
    b = (double *)aligned_malloc(sizeof(double)*sz_trid[0]*sz_trid[1]*sz_trid[2] + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
    double *b_tmp = (double *)malloc(sz_trid[2] * sizeof(double));
    double *dzf = (double *)malloc((sz_trid[2] + 2) * sizeof(double));
    
    for (int i = 0; i < sz_trid[2] + 2; i++) {
        dzf[i] = dzf_global[st_trid[2]-1+i];
    }
    set_trid_coeff(sz_trid[2], dzf, bctype_z, 'c', &neighbor_trid[4], a, b_tmp, c);
    #pragma omp parallel
    {
        for (int k = 0; k < sz_trid[2]; k++) {
            #pragma omp for schedule(static)
            for (int j = 0; j < sz_trid[1]; j++) {
                for (int i = 0; i < sz_trid[0]; i++) {
                    b[i + j * sz_trid[0] + k * sz_trid[0] * sz_trid[1]] = 
                        b_tmp[k] + lambdaxy[i + j * sz_trid[0]];
                }
            }
        }
        // decompose coefficient b
        for (int k = 1; k < sz_trid[2]; k++) {
            #pragma omp for schedule(static)
            for (int j = 0; j < sz_trid[1]; j++) {
                for (int i = 0; i < sz_trid[0]; i++) {
                    double a_tmp = a[k] / b[i + j * sz_trid[0] + (k-1) * sz_trid[0] * sz_trid[1]];
                    b[i + j * sz_trid[0] + k * sz_trid[0] * sz_trid[1]] -= (a_tmp * c[k-1]);
                }
            }
        }
    }
    // double *b_k_ptr = b + sz_trid[0] * sz_trid[1];
    // double *b_km1_ptr = b;
    // for (int k = 1; k < sz_trid[2]; k++) {
    //     for (int j = 0; j < sz_trid[1]; j++) {
    //         for (int i = 0; i < sz_trid[0]; i++) {
    //             double a_tmp = a[k] / *b_km1_ptr;
    //             *b_k_ptr -= a_tmp * c[k-1];
    //             b_k_ptr++;
    //             b_km1_ptr++;
    //         }
    //     }
    // }

    // determine whether the tridiagonal systems are periodic or not
    // NOTE: not yet implemented

    // calculate the correction of the right-hand side according to BC types in x, y, z direction
    // NOTE: not yet implemented

#ifdef _PDD
    COMM_CART_PDD = decomp_2d_comm_cart_y;
    
    // initialize work arrays for PDD algorithm
    init_pdd_array(sz_trid);
#endif

    free(lambdax);
    free(lambday);
    free(lambdaxy);
    free(b_tmp);
    free(dzf);

}

void ssolve_trid(int *sz, double *var) {

    int ie = sz[0];
    int je = sz[1];
    int ke = sz[2];

    #pragma omp parallel
    {
        for (int k = 1; k < ke; k++) {
            #pragma omp for schedule(static)
            for (int j = 0; j < je; j++) {
                double *var_down = &var[k*sz[0]*sz[1]+j*sz[0]];
                double *var_up = &var[(k-1)*sz[0]*sz[1]+j*sz[0]];
                double *b_p = &b[(k-1)*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
                for (int i = 0; i < ie; i++) {
                    double a_tmp = a[k]/b_p[i];
                    var_down[i] -= a_tmp * var_up[i];
                }
            }
        }

        #pragma omp for schedule(static)
        for (int j = 0; j < je; j++) {
            double *b_p = &b[(ke-1)*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
            double *var_p = &var[(ke-1)*sz[0]*sz[1]+j*sz[0]];
            for (int i = 0; i < ie; i++) {
                if (b_p[i] != 0) {
                    var_p[i] /= b_p[i];
                }
                else {
                    var_p[i] = 0;
                }
            }
        }
        
        for (int k = ke-2; k >=0 ; k--) {
            #pragma omp for schedule(static)
            for (int j = 0; j < je; j++) {
                double *var_down = &var[k*sz[0]*sz[1]+j*sz[0]];
                double *var_up = &var[(k+1)*sz[0]*sz[1]+j*sz[0]];
                double *b_p = &b[k*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
                for (int i = 0; i < ie; i++) {
                    var_down[i] = (var_down[i]-c[k]*var_up[i])/b_p[i];
                }
            }
        }
    }
    
}

#ifdef _PDD
void psolve_trid(int *sz, double *var) {

    int ie = sz[0];
    int je = sz[1];
    int ke = sz[2];

    #pragma omp parallel
    {
        for (int k = 1; k < ke; k++) {
            #pragma omp for schedule(static)
            for (int j = 0; j < je; j++) {
                double *var_down = &var[k*sz[0]*sz[1]+j*sz[0]];
                double *var_up = &var[(k-1)*sz[0]*sz[1]+j*sz[0]];
                double *b_p = &b[(k-1)*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
                for (int i = 0; i < ie; i++) {
                    double a_tmp = a[k]/b_p[i];
                    var_down[i] -= a_tmp * var_up[i];
                }
            }
        }

        #pragma omp for schedule(static)
        for (int j = 0; j < je; j++) {
            double *b_p = &b[(ke-1)*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
            double *var_p = &var[(ke-1)*sz[0]*sz[1]+j*sz[0]];
            for (int i = 0; i < ie; i++) {
                if (b_p[i] != 0) {
                    var_p[i] /= b_p[i];
                }
                else {
                    var_p[i] = 0;
                }
            }
        }
        
        for (int k = ke-2; k >=0 ; k--) {
            #pragma omp for schedule(static)
            for (int j = 0; j < je; j++) {
                double *var_up = &var[k*sz[0]*sz[1]+j*sz[0]];
                double *var_down = &var[(k+1)*sz[0]*sz[1]+j*sz[0]];
                double *b_p = &b[k*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
                for (int i = 0; i < ie; i++) {
                    var_up[i] = (var_up[i]-c[k]*var_down[i])/b_p[i];
                }
            }
        }
    }
    

    #ifdef GPTL
        GPTLstart("----Comm in PDD");
    #endif

    MPI_Sendrecv(var, tmp_var_pdd_len, MPI_DOUBLE, neighbor_trid[4], 100,
         tmp_var_pdd, tmp_var_pdd_len, MPI_DOUBLE, neighbor_trid[5], 100,
         COMM_CART_PDD, MPI_STATUS_IGNORE);

    #ifdef GPTL
        GPTLstop("----Comm in PDD");
    #endif

    if (neighbor_trid[5] != MPI_PROC_NULL) {
        #pragma omp parallel for schedule(static)
        for (int j = 0; j < je; j++) {
            double *w_pdd_p = &w_pdd[(ke-1)*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
            double *tmp_v_pdd_p = &tmp_v_pdd[j*sz_trid[0]];
            double *var_p = &var[(ke-1)*sz[0]*sz[1]+j*sz[0]];
            double *y2_pdd_p = &y2_pdd[j*sz_trid[0]];
            double *y3_pdd_p = &y3_pdd[j*sz_trid[0]];
            double *tmp_var_pdd_p = &tmp_var_pdd[j*sz_trid[0]];
            for (int i = 0; i < ie; i++) {
                double det_pdd = w_pdd_p[i] * tmp_v_pdd_p[i] - ((double)1.0);
                y2_pdd_p[i] = (var_p[i] * tmp_v_pdd_p[i] - tmp_var_pdd_p[i]) / det_pdd;
                y3_pdd_p[i] = (tmp_var_pdd_p[i] * w_pdd_p[i] - var_p[i]) / det_pdd;
            }
        }
    }

    #ifdef GPTL
        GPTLstart("----Comm in PDD");
    #endif

    MPI_Sendrecv(y3_pdd, tmp_var_pdd_len, MPI_DOUBLE, neighbor_trid[5], 100,
                 y1_pdd, tmp_var_pdd_len, MPI_DOUBLE, neighbor_trid[4], 100,
                 COMM_CART_PDD, MPI_STATUS_IGNORE);

    #ifdef GPTL
        GPTLstop("----Comm in PDD");
    #endif

    #pragma omp parallel for schedule(static)
    for (int k = 0; k < ke; k++) {
        for (int j = 0; j < je; j++) {
            double *var_p = &var[k*sz[0]*sz[1]+j*sz[0]];
            double *v_pdd_p = &v_pdd[k*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
            double *w_pdd_p = &w_pdd[k*sz_trid[0]*sz_trid[1]+j*sz_trid[0]];
            double *y1_pdd_p = &y1_pdd[j*sz_trid[0]];
            double *y2_pdd_p = &y2_pdd[j*sz_trid[0]];
            for (int i = 0; i < ie; i++) {
                var_p[i] = var_p[i] - v_pdd_p[i]*y1_pdd_p[i] - w_pdd_p[i]*y2_pdd_p[i];
            }
        }
    }
}
#endif

void execute_poisson_solver(double *var) {

#ifdef USE_MKL
    static fft_mkl_plan_ptr plan_fwd_xpen, plan_bwd_xpen, plan_fwd_ypen, plan_bwd_ypen;
#else
    static fftw_plan plan_fwd_xpen, plan_bwd_xpen, plan_fwd_ypen, plan_bwd_ypen;
#endif
    plan_fwd_xpen = plan[0][0];
    plan_bwd_xpen = plan[0][1];
    plan_fwd_ypen = plan[1][0];
    plan_bwd_ypen = plan[1][1];

    #ifdef GPTL
        GPTLstart("--Copy & Fwd X-FFT");
    #endif

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < xsize[2]; k++) {
        double *dst = &var_xpen[k*var_xpen_offset1];
        double *src = &var[(k+1)*var_offset1+var_offset2+1];
        for (size_t j = 0; j < xsize[1]; j++) {
            memcpy(dst, src, xsize[0] * sizeof(double));
            dst += var_xpen_offset2;
            src += var_offset2;
        }
        execute_fft(plan_fwd_xpen, &var_xpen[k*var_xpen_offset1]);
    }

    #ifdef GPTL
        GPTLstop("--Copy & Fwd X-FFT");
        GPTLstart("--Transpose x to y");
    #endif

    if (var_xpen != var_ypen) {
        transpose_x_to_y_real(var_xpen, var_ypen, xsize, ysize, &decomp_main);
    }

    #ifdef GPTL
        GPTLstop("--Transpose x to y");
        GPTLstart("--Forward Y-FFT");
    #endif

    size_t var_ypen_offset1 = ysize[0] * ysize[1];
    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < ysize[2]; k++) {
        execute_fft(plan_fwd_ypen, &var_ypen[k*var_ypen_offset1]);
    }

    #ifdef GPTL
        GPTLstop("--Forward Y-FFT");
    #endif

    #ifdef _PDD
    {
        #ifdef GPTL
            GPTLstart("--Solve trid");
        #endif

        psolve_trid(ysize, var_ypen);

        #ifdef GPTL
            GPTLstop("--Solve trid");
        #endif
    }
    #else
    {
        #ifdef GPTL
            GPTLstart("--Transpose y to z");
        #endif

        if (var_ypen != var_zpen) {
            transpose_y_to_z_real(var_ypen, var_zpen, ysize, zsize, &decomp_main);
        }

        #ifdef GPTL
            GPTLstop("--Transpose y to z");
            GPTLstart("--Solve trid");
        #endif

        ssolve_trid(zsize, var_zpen);

        #ifdef GPTL
            GPTLstop("--Solve trid");
            GPTLstart("--Transpose z to y");
        #endif

        if (var_ypen != var_zpen) {
            transpose_z_to_y_real(var_zpen, var_ypen, zsize, ysize, &decomp_main);
        }

        #ifdef GPTL
            GPTLstop("--Transpose z to y");
        #endif
    }
    #endif

    #ifdef GPTL
        GPTLstart("--Backward Y-FFT");
    #endif

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < ysize[2]; k++) {
        execute_fft(plan_bwd_ypen, &var_ypen[k*var_ypen_offset1]);
    }

    #ifdef GPTL
        GPTLstop("--Backward Y-FFT");
        GPTLstart("--Transpose y to x");
    #endif

    if (var_xpen != var_ypen) {
        transpose_y_to_x_real(var_ypen, var_xpen, ysize, xsize, &decomp_main);
    }

    #ifdef GPTL
        GPTLstop("--Transpose y to x");
        GPTLstart("--Bwd X-FFT & Copy");
    #endif

    #pragma omp parallel for schedule(static)
    for (size_t k = 0; k < xsize[2]; k++) {
        double *src = &var_xpen[k*var_xpen_offset1];
        double *dst = &var[(k+1)*var_offset1+var_offset2+1];
        execute_fft(plan_bwd_xpen, src);
        for (size_t j = 0; j < xsize[1]; j++) {
            for (size_t i = 0; i < xsize[0]; i++) {
                dst[i] = src[i] * fft_normfactor;
            }
            src += var_xpen_offset2;
            dst += var_offset2;
        }
    }

    #ifdef GPTL
        GPTLstop("--Bwd X-FFT & Copy");
    #endif

}

void free_poisson_solver() {
    
    // release work arrays
    if (var_xpen != NULL) aligned_free(var_xpen);
    if (var_xpen != var_ypen && var_ypen != NULL) aligned_free(var_ypen);
    if (var_ypen != var_zpen && var_zpen != NULL) aligned_free(var_zpen);

    // release tridiagonal coefficients arrays
    if (a != NULL) aligned_free(a);
    if (b != NULL) aligned_free(b);
    if (c != NULL) aligned_free(c);

#ifdef _PDD
    // release PDD related arrays
    if (v_pdd != NULL) aligned_free(v_pdd);
    if (w_pdd != NULL) aligned_free(w_pdd);
    if (y1_pdd != NULL) aligned_free(y1_pdd);
    if (y2_pdd != NULL) aligned_free(y2_pdd);
    if (y3_pdd != NULL) aligned_free(y3_pdd);
    if (tmp_v_pdd != NULL) aligned_free(tmp_v_pdd);
    if (tmp_var_pdd != NULL) aligned_free(tmp_var_pdd);
#endif

    // release FFT
    free_fft(plan);
}
