#ifndef DECOMP_2D_H
#define DECOMP_2D_H

#include <stdbool.h>
#include <mpi.h>

#ifdef SINGLE_PRECISION
typedef float DECOMP_2D_REAL;
#define DECOMP_2D_MPI_REAL MPI_FLOAT
#else
typedef double DECOMP_2D_REAL;
#define DECOMP_2D_MPI_REAL MPI_DOUBLE
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern int decomp_2d_nx_global, decomp_2d_ny_global, decomp_2d_nz_global;
extern int decomp_2d_periodic_x, decomp_2d_periodic_y, decomp_2d_periodic_z;
extern int decomp_2d_nrank, decomp_2d_nproc;
extern MPI_Comm decomp_2d_comm_cart_x;
extern MPI_Comm decomp_2d_comm_cart_y;
extern MPI_Comm decomp_2d_comm_cart_z;
extern MPI_Fint decomp_2d_comm_cart_x_f;
extern MPI_Fint decomp_2d_comm_cart_y_f;
extern MPI_Fint decomp_2d_comm_cart_z_f;
extern MPI_Comm decomp_2d_comm_row;
extern MPI_Comm decomp_2d_comm_col;

typedef struct DECOMP_2D_INFO {
    // send/receive buffer counts and displacements for MPI_Alltoallv
    int *x1cnts, *x1disp;
    int *y1cnts, *y1disp;
    int *y2cnts, *y2disp;
    int *z2cnts, *z2disp;
    // send/receive buffer counts for MPI_Alltoall: either for evenly distributed data or for padded-Alltoall
    int x1count, y1count, y2count, z2count;

    // 2D decomposition information
    // staring/ending index and size of data held by current rank
    int xst[3], xen[3], xsz[3]; // x-pencil
    int yst[3], yen[3], ysz[3]; // y-pencil
    int zst[3], zen[3], zsz[3]; // z-pencil
    
    // how each dimension is distributed along pencils
    int *x1st, *x1en, *x1dist;  // distribute nx_global grid points to dims[0] ranks
    int *y1st, *y1en, *y1dist;  // distribute ny_global grid points to dims[0] ranks
    int *y2st, *y2en, *y2dist;  // distribute ny_global grid points to dims[1] ranks
    int *z2st, *z2en, *z2dist;  // distribute nz_global grid points to dims[1] ranks

    // evenly distributed data
    bool even;
} decomp_2d_info;
extern decomp_2d_info decomp_main;

extern int dims[2];

extern void decomp_2d_init(int nx, int ny, int nz, int p_row, int p_col, bool periodic_bc[3]);
extern void decomp_2d_finalize();
extern void decomp_2d_abort(int errorcode, char *msg);
extern void decomp_info_init(decomp_2d_info *decomp, int nx, int ny, int nz);
extern void decomp_info_finalize(decomp_2d_info *decomp);
extern size_t get_decomp_2d_work_size();
extern DECOMP_2D_REAL *get_decomp_2d_work1_r();
extern DECOMP_2D_REAL *get_decomp_2d_work2_r();
extern void transpose_x_to_y_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp);
extern void transpose_y_to_x_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp);
extern void transpose_y_to_z_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp);
extern void transpose_z_to_y_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp);

#ifdef __cplusplus
}
#endif

#endif // DECOMP_2D_H
