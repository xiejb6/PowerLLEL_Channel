#include "memory.h"
#include "decomp_2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#ifdef GPTL
#include <gptl.h>
#endif

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

// public global variables definitions
int decomp_2d_nx_global, decomp_2d_ny_global, decomp_2d_nz_global;
int decomp_2d_periodic_x, decomp_2d_periodic_y, decomp_2d_periodic_z;
int decomp_2d_nrank, decomp_2d_nproc;
MPI_Comm decomp_2d_comm_cart_x;
MPI_Comm decomp_2d_comm_cart_y;
MPI_Comm decomp_2d_comm_cart_z;
MPI_Fint decomp_2d_comm_cart_x_f;
MPI_Fint decomp_2d_comm_cart_y_f;
MPI_Fint decomp_2d_comm_cart_z_f;
MPI_Comm decomp_2d_comm_row;
MPI_Comm decomp_2d_comm_col;
decomp_2d_info decomp_main;

// private global variables - parameters for MPI 2D Cartesian topology
int dims[2];
static int coords[2];

// private global variables - internal buffers used by MPI_Alltoall(v)
static size_t decomp_buf_size = 0;
static DECOMP_2D_REAL *work1_r = NULL;
static DECOMP_2D_REAL *work2_r = NULL;

static void partition(int n[3], int pdim[3], int coords[2], 
    int dim1st[], int dim1en[], int dim1sz[], 
    int dim2st[], int dim2en[], int dim2sz[], 
    int lstart[3], int lend[3], int lsize[3]);
static void prepare_buffer(decomp_2d_info *decomp);
static void distribute(int num, int nprocs, int st[], int en[], int sz[]);
static void distribute_among_threads(char c_or_f, int nthreads, int tid, int num, int factor, int *sz, int *st, int *en);

void decomp_2d_init(int nx, int ny, int nz, int p_row, int p_col, bool periodic_bc[3])
{
    decomp_2d_nx_global = nx;
    decomp_2d_ny_global = ny;
    decomp_2d_nz_global = nz;
    decomp_2d_periodic_x = periodic_bc[0] ? 1 : 0;
    decomp_2d_periodic_y = periodic_bc[1] ? 1 : 0;
    decomp_2d_periodic_z = periodic_bc[2] ? 1 : 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &decomp_2d_nrank);
    MPI_Comm_size(MPI_COMM_WORLD, &decomp_2d_nproc);

    if (decomp_2d_nproc != p_row*p_col) {
        decomp_2d_abort(1, "Invalid 2D processor grid - nproc /= p_row*p_col");
    } else {
        dims[0] = p_row;
        dims[1] = p_col;
    }

    // Create 2D Catersian topology
    // Note that in order to support periodic B.C. in the halo-cell code,
    // need to create multiple topology objects: decomp_2d_comm_cart_?,
    // corresponding to three pencil orientations. They contain almost
    // identical topological information but allow different combinations
    // of periodic conditions.
    int periods[2];
    periods[0] = decomp_2d_periodic_y;
    periods[1] = decomp_2d_periodic_z;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &decomp_2d_comm_cart_x);
    periods[0] = decomp_2d_periodic_x;
    periods[1] = decomp_2d_periodic_z;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &decomp_2d_comm_cart_y);
    periods[0] = decomp_2d_periodic_x;
    periods[1] = decomp_2d_periodic_y;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &decomp_2d_comm_cart_z);

    decomp_2d_comm_cart_x_f = MPI_Comm_c2f(decomp_2d_comm_cart_x);
    decomp_2d_comm_cart_y_f = MPI_Comm_c2f(decomp_2d_comm_cart_y);
    decomp_2d_comm_cart_z_f = MPI_Comm_c2f(decomp_2d_comm_cart_z);

    MPI_Cart_coords(decomp_2d_comm_cart_x, decomp_2d_nrank, 2, coords);

    // derive communicators defining sub-groups for Alltoall(v)
    int remain_dims[2] = {1, 0};
    MPI_Cart_sub(decomp_2d_comm_cart_x, remain_dims, &decomp_2d_comm_col);
    remain_dims[0] = 0;
    remain_dims[1] = 1;
    MPI_Cart_sub(decomp_2d_comm_cart_x, remain_dims, &decomp_2d_comm_row);

    // actually generate all 2D decomposition information
    decomp_info_init(&decomp_main, nx, ny, nz);

#ifdef EVEN
    if (decomp_2d_nrank == 0)
        printf("Padded ALLTOALL optimisation on\n");
#endif
}

void decomp_2d_finalize()
{
    MPI_Comm_free(&decomp_2d_comm_cart_x);
    MPI_Comm_free(&decomp_2d_comm_cart_y);
    MPI_Comm_free(&decomp_2d_comm_cart_z);
    MPI_Comm_free(&decomp_2d_comm_row);
    MPI_Comm_free(&decomp_2d_comm_col);
    
    decomp_info_finalize(&decomp_main);
}

void decomp_2d_abort(int errorcode, char *msg)
{
    int nrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
    if (nrank == 0)
        printf("DECOMP_2D ERROR - errorcode: %d, message: %s \n", errorcode, msg);
    MPI_Abort(MPI_COMM_WORLD, errorcode);
}

void decomp_info_init(decomp_2d_info *decomp, int nx, int ny, int nz)
{
    // verify the global size can actually be distributed as pencils
    if (nx < dims[0] || 
        ny < dims[0] || 
        ny < dims[1] || 
        nz < dims[1]) {
        decomp_2d_abort(6, 
            "Invalid 2D processor grid. Make sure that min(nx,ny) >= p_row and min(ny,nz) >= p_col");
    }

    if (nx % dims[0] == 0 &&
        ny % dims[0] == 0 &&
        ny % dims[1] == 0 &&
        nz % dims[1] == 0) {
        decomp->even = true;
    } else {
        decomp->even = false;
    }

    // distribute grid points across processors
    //   e.g. 17 points across 4 processors would be distibuted as (4,4,4,5)
    //   such global information is required locally at MPI_Alltoall(v) time
    decomp->x1st   = (int *) malloc(dims[0] * sizeof(int));
    decomp->x1en   = (int *) malloc(dims[0] * sizeof(int));
    decomp->x1dist = (int *) malloc(dims[0] * sizeof(int));
    distribute(nx, dims[0], decomp->x1st, decomp->x1en, decomp->x1dist);
    decomp->y1st   = (int *) malloc(dims[0] * sizeof(int));
    decomp->y1en   = (int *) malloc(dims[0] * sizeof(int));
    decomp->y1dist = (int *) malloc(dims[0] * sizeof(int));
    distribute(ny, dims[0], decomp->y1st, decomp->y1en, decomp->y1dist);
    decomp->y2st   = (int *) malloc(dims[1] * sizeof(int));
    decomp->y2en   = (int *) malloc(dims[1] * sizeof(int));
    decomp->y2dist = (int *) malloc(dims[1] * sizeof(int));
    distribute(ny, dims[1], decomp->y2st, decomp->y2en, decomp->y2dist);
    decomp->z2st   = (int *) malloc(dims[1] * sizeof(int));
    decomp->z2en   = (int *) malloc(dims[1] * sizeof(int));
    decomp->z2dist = (int *) malloc(dims[1] * sizeof(int));
    distribute(nz, dims[1], decomp->z2st, decomp->z2en, decomp->z2dist);

    // generate partition information - staring/ending index and size of data held by current rank
    int n[3] = {nx, ny, nz};
    int pdim[3];
    pdim[0] = 0; pdim[1] = 1; pdim[2] = 2;
    partition(n, pdim, coords, 
              decomp->y1st, decomp->y1en, decomp->y1dist,
              decomp->z2st, decomp->z2en, decomp->z2dist,
              decomp->xst, decomp->xen, decomp->xsz);
    pdim[0] = 1; pdim[1] = 0; pdim[2] = 2;
    partition(n, pdim, coords, 
              decomp->x1st, decomp->x1en, decomp->x1dist,
              decomp->z2st, decomp->z2en, decomp->z2dist,
              decomp->yst, decomp->yen, decomp->ysz);
    pdim[0] = 1; pdim[1] = 2; pdim[2] = 0;
    partition(n, pdim, coords, 
              decomp->x1st, decomp->x1en, decomp->x1dist,
              decomp->y2st, decomp->y2en, decomp->y2dist,
              decomp->zst, decomp->zen, decomp->zsz);
    
    // prepare send/receive buffer displacement and count for Alltoall(v)
    decomp->x1cnts = (int *) malloc(dims[0] * sizeof(int));
    decomp->x1disp = (int *) malloc(dims[0] * sizeof(int));
    decomp->y1cnts = (int *) malloc(dims[0] * sizeof(int));
    decomp->y1disp = (int *) malloc(dims[0] * sizeof(int));
    decomp->y2cnts = (int *) malloc(dims[1] * sizeof(int));
    decomp->y2disp = (int *) malloc(dims[1] * sizeof(int));
    decomp->z2cnts = (int *) malloc(dims[1] * sizeof(int));
    decomp->z2disp = (int *) malloc(dims[1] * sizeof(int));
    prepare_buffer(decomp);

    // allocate memory for the MPI_Alltoall(v) buffers
    size_t buf_size = decomp->xsz[0] * decomp->xsz[1] * decomp->xsz[2];
    buf_size = MAX( buf_size, decomp->ysz[0] * decomp->ysz[1] * decomp->ysz[2] );
    buf_size = MAX( buf_size, decomp->zsz[0] * decomp->zsz[1] * decomp->zsz[2] );
#ifdef EVEN
    buf_size = MAX( buf_size, decomp->x1count * dims[0] );
    buf_size = MAX( buf_size, decomp->y2count * dims[1] );
#endif
    if (buf_size > decomp_buf_size) {
        decomp_buf_size = buf_size;
        if (work1_r != NULL) aligned_free(work1_r);
        if (work2_r != NULL) aligned_free(work2_r);
        work1_r = (DECOMP_2D_REAL *) aligned_malloc(buf_size * sizeof(DECOMP_2D_REAL) + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
        work2_r = (DECOMP_2D_REAL *) aligned_malloc(buf_size * sizeof(DECOMP_2D_REAL) + MEM_ALIGN_SIZE, MEM_ALIGN_SIZE);
        if (work1_r == NULL || work2_r == NULL)
            decomp_2d_abort(2, "Out of memory when allocating 2DECOMP workspace");
    }
    
}

void decomp_info_finalize(decomp_2d_info *decomp)
{
    if (decomp->x1st   != NULL) free(decomp->x1st  );
    if (decomp->x1en   != NULL) free(decomp->x1en  );
    if (decomp->x1dist != NULL) free(decomp->x1dist);
    if (decomp->y1st   != NULL) free(decomp->y1st  );
    if (decomp->y1en   != NULL) free(decomp->y1en  );
    if (decomp->y1dist != NULL) free(decomp->y1dist);
    if (decomp->y2st   != NULL) free(decomp->y2st  );
    if (decomp->y2en   != NULL) free(decomp->y2en  );
    if (decomp->y2dist != NULL) free(decomp->y2dist);
    if (decomp->z2st   != NULL) free(decomp->z2st  );
    if (decomp->z2en   != NULL) free(decomp->z2en  );
    if (decomp->z2dist != NULL) free(decomp->z2dist);

    if (decomp->x1cnts != NULL) free(decomp->x1cnts);
    if (decomp->x1disp != NULL) free(decomp->x1disp);
    if (decomp->y1cnts != NULL) free(decomp->y1cnts);
    if (decomp->y1disp != NULL) free(decomp->y1disp);
    if (decomp->y2cnts != NULL) free(decomp->y2cnts);
    if (decomp->y2disp != NULL) free(decomp->y2disp);
    if (decomp->z2cnts != NULL) free(decomp->z2cnts);
    if (decomp->z2disp != NULL) free(decomp->z2disp);

    decomp_buf_size = 0;
    if (work1_r != NULL) aligned_free(work1_r);
    if (work2_r != NULL) aligned_free(work2_r);
}

size_t get_decomp_2d_work_size() {
    return decomp_buf_size * sizeof(DECOMP_2D_REAL);
}

DECOMP_2D_REAL *get_decomp_2d_work1_r() {
    return work1_r;
}

DECOMP_2D_REAL *get_decomp_2d_work2_r() {
    return work2_r;
}

// Distributes grid points in one dimension, handles uneven distribution properly
// Note the last blocks along pencils always get assigned more grid points
//   num    -- data size in any dimension to be partitioned
//   nprocs -- number of processors in that dimension
//   st     -- array of starting index
//   en     -- array of ending index
//   sz     -- array of local size  (redundent)
void distribute(int num, int nprocs, int st[], int en[], int sz[])
{
    int size = num / nprocs;
    int nu = num - size * nprocs;
    int nl = nprocs - nu;
    st[0] = 1;
    sz[0] = size;
    en[0] = size;
    for (int i = 1; i < nl; i++) {
        st[i] = st[i-1] + size;
        sz[i] = size;
        en[i] = en[i-1] + size;
    }
    size++;
    for (int i = nl; i < nprocs; i++)
    {
        st[i] = en[i-1] + 1;
        sz[i] = size;
        en[i] = en[i-1] + size;
    }
    en[nprocs-1] = num;
    sz[nprocs-1] = num - st[nprocs-1] + 1;
}

// Distributes grid points among threads
// Note that `c_or_f` determines value of `st`, in C/Fortran style
//   e.g., distribute 7 among 3 threads (factor = 1)
//   thread :       0           1          2
//   c      : [0, 1, 2, 3), [3, 4, 5), [5, 6, 7)
//   fortran: [1, 2, 3],    [4, 5],    [6, 7] 
void distribute_among_threads(char c_or_f, int nthreads, int tid, int num, int factor, int *sz, int *st, int *en)
{
    *sz = num / nthreads;
    int nloops = *sz/factor;
    *sz = nloops * factor;
    int remainder = num - (*sz) * nthreads;
    if (tid < remainder/factor) {
        *sz += factor;
        *st = tid * (*sz);
    } else {
        *st = remainder + tid * (*sz);
    }
    *en = *st + *sz;

    if (c_or_f == 'f') *st += 1;
}

// Find sub-domain information held by current processor
//   pdim[3] - number of processor grid in each dimension, valid values:
//     0 - distribute locally
//     1 - distribute across p_row(or dims[0])
//     2 - distribute across p_col(or dims[1])
static void partition(int n[3], int pdim[3], int coords[2], 
    int dim1st[], int dim1en[], int dim1sz[], 
    int dim2st[], int dim2en[], int dim2sz[], 
    int lstart[3], int lend[3], int lsize[3])
{
    for (int i = 0; i < 3; i++) {
        int gsize = n[i];
        switch (pdim[i])
        {
        case 0:
            lstart[i] = 1;
            lend  [i] = gsize;
            lsize [i] = gsize;
            break;
        case 1:
            lstart[i] = dim1st[coords[0]];
            lend  [i] = dim1en[coords[0]];
            lsize [i] = dim1sz[coords[0]];
            break;
        case 2:
            lstart[i] = dim2st[coords[1]];
            lend  [i] = dim2en[coords[1]];
            lsize [i] = dim2sz[coords[1]];
            break;
        }
    }
}

// Prepare the send / receive buffers for MPI_Alltoall(v) communications
static void prepare_buffer(decomp_2d_info *decomp)
{
    // MPI_Alltoallv buffer information
    for (int i = 0; i < dims[0]; i++) {
        decomp->x1cnts[i] = decomp->x1dist[i] * decomp->xsz[1] * decomp->xsz[2];
        decomp->y1cnts[i] = decomp->ysz[0] * decomp->y1dist[i] * decomp->ysz[2];
        if (i == 0) {
            decomp->x1disp[i] = 0;
            decomp->y1disp[i] = 0;
        } else {
            decomp->x1disp[i] = decomp->x1disp[i-1] + decomp->x1cnts[i-1];
            decomp->y1disp[i] = decomp->y1disp[i-1] + decomp->y1cnts[i-1];
        }
    }

    for (int i = 0; i < dims[1]; i++) {
        decomp->y2cnts[i] = decomp->ysz[0] * decomp->y2dist[i] * decomp->ysz[2];
        decomp->z2cnts[i] = decomp->zsz[0] * decomp->zsz[1] * decomp->z2dist[i];
        if (i == 0) {
            decomp->y2disp[i] = 0;
            decomp->z2disp[i] = 0;
        } else {
            decomp->y2disp[i] = decomp->y2disp[i-1] + decomp->y2cnts[i-1];
            decomp->z2disp[i] = decomp->z2disp[i-1] + decomp->z2cnts[i-1];
        }
    }

    // MPI_Alltoall buffer information
    //   For unevenly distributed data, pad smaller messages. Note the 
    //   last blocks along pencils always get assigned more grid points.
    //   For evenly distributed data, no padding.
    decomp->x1count = decomp->x1dist[ dims[0]-1 ] *
                      decomp->y1dist[ dims[0]-1 ] *
                      decomp->xsz[2];
    decomp->y1count = decomp->x1count;
    decomp->y2count = decomp->ysz[0] *
                      decomp->y2dist[ dims[1]-1 ] *
                      decomp->z2dist[ dims[1]-1 ];
    decomp->z2count = decomp->y2count;
}


// ===============================
// Pencil transposition functions 
// ===============================

// Pack Alltoall(v) buffers for x_to_y transpose
static void mem_split_xy_real(DECOMP_2D_REAL *in, int n1, int n2, int n3, DECOMP_2D_REAL *out, int iproc, int *x1st, int *x1en, int *x1disp) {
    int in_offset1 = n1 * n2;
    int in_offset2 = n1;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int m = 0; m < iproc; m++) {
        int i1 = x1st[m]-1;
        int i2 = x1en[m];
        int pos = x1disp[m];

        for (int k = 0; k < n3; k++) {
            for (int j = 0; j < n2; j++) {
                for (int i = i1; i < i2; i++) {
                    out[pos] = in[i + j * in_offset2 + k * in_offset1];
                    pos++;
                }
            }
        }
    }
}

// Pack Alltoall(v) buffers for y_to_x transpose
static void mem_split_yx_real(DECOMP_2D_REAL *in, int n1, int n2, int n3, DECOMP_2D_REAL *out, int iproc, int *y1st, int *y1en, int *y1disp) {
    int in_offset1 = n1 * n2;
    int in_offset2 = n1;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int m = 0; m < iproc; m++) {
        int i1 = y1st[m]-1;
        int i2 = y1en[m];
        int pos = y1disp[m];

        for (int k = 0; k < n3; k++) {
            for (int j = i1; j < i2; j++) {
                for (int i = 0; i < n1; i++) {
                    out[pos] = in[i + j * in_offset2 + k * in_offset1];
                    pos++;
                }
            }
        }
    }
}

// Pack Alltoall(v) buffers for y_to_z transpose
static void mem_split_yz_real(DECOMP_2D_REAL *in, int n1, int n2, int n3, DECOMP_2D_REAL *out, int iproc, int *y2st, int *y2en, int *y2disp) {
    int in_offset1 = n1 * n2;
    int in_offset2 = n1;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int m = 0; m < iproc; m++) {
        int i1 = y2st[m]-1;
        int i2 = y2en[m];
        int pos = y2disp[m];

        for (int k = 0; k < n3; k++) {
            for (int j = i1; j < i2; j++) {
                for (int i = 0; i < n1; i++) {
                    out[pos] = in[i + j * in_offset2 + k * in_offset1];
                    pos++;
                }
            }
        }
    }
}

// Unpack Alltoall(v) buffers for x_to_y transpose
static void mem_merge_xy_real(DECOMP_2D_REAL *in, int n1, int n2, int n3, DECOMP_2D_REAL *out, int iproc, int *y1st, int *y1en, int *y1disp) {
    int out_offset1 = n1 * n2;
    int out_offset2 = n1;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int m = 0; m < iproc; m++) {
        int i1 = y1st[m]-1;
        int i2 = y1en[m];
        int pos = y1disp[m];

        for (int k = 0; k < n3; k++) {
            for (int j = i1; j < i2; j++) {
                for (int i = 0; i < n1; i++) {
                    out[i + j * out_offset2 + k * out_offset1] = in[pos];
                    pos++;
                }
            }
        }
    }
}

// Unpack Alltoall(v) buffers for y_to_x transpose
static void mem_merge_yx_real(DECOMP_2D_REAL *in, int n1, int n2, int n3, DECOMP_2D_REAL *out, int iproc, int *x1st, int *x1en, int *x1disp) {
    int out_offset1 = n1 * n2;
    int out_offset2 = n1;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int m = 0; m < iproc; m++) {
        int i1 = x1st[m]-1;
        int i2 = x1en[m];
        int pos = x1disp[m];

        for (int k = 0; k < n3; k++) {
            for (int j = 0; j < n2; j++) {
                for (int i = i1; i < i2; i++) {
                    out[i + j * out_offset2 + k * out_offset1] = in[pos];
                    pos++;
                }
            }
        }
    }
}

// Unpack Alltoall(v) buffers for z_to_y transpose
static void mem_merge_zy_real(DECOMP_2D_REAL *in, int n1, int n2, int n3, DECOMP_2D_REAL *out, int iproc, int *y2st, int *y2en, int *y2disp) {
    int out_offset1 = n1 * n2;
    int out_offset2 = n1;

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int m = 0; m < iproc; m++) {
        int i1 = y2st[m]-1;
        int i2 = y2en[m];
        int pos = y2disp[m];

        for (int k = 0; k < n3; k++) {
            for (int j = i1; j < i2; j++) {
                for (int i = 0; i < n1; i++) {
                    out[i + j * out_offset2 + k * out_offset1] = in[pos];
                    pos++;
                }
            }
        }
    }
}

// Transpose data (real) from X to Y pencil
void transpose_x_to_y_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp) {

    #ifdef GPTL
        GPTLstart("----mem_split_xy");
    #endif
    
    mem_split_xy_real(src, src_size[0], src_size[1], src_size[2], work1_r, dims[0], decomp->x1st, decomp->x1en, decomp->x1disp);

    #ifdef GPTL
        GPTLstop("----mem_split_xy");
        GPTLstart("----mpi_alltoall_xy");
    #endif

    MPI_Alltoallv(work1_r, decomp->x1cnts, decomp->x1disp, DECOMP_2D_MPI_REAL, 
                  work2_r, decomp->y1cnts, decomp->y1disp, DECOMP_2D_MPI_REAL,
                  decomp_2d_comm_col);

    #ifdef GPTL
        GPTLstop("----mpi_alltoall_xy");
        GPTLstart("----mem_merge_xy");
    #endif

    mem_merge_xy_real(work2_r, dst_size[0], dst_size[1], dst_size[2], dst, dims[0], decomp->y1st, decomp->y1en, decomp->y1disp);

    #ifdef GPTL
        GPTLstop("----mem_merge_xy");
    #endif
}

// Transpose data (real) from Y to X pencil
void transpose_y_to_x_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp) {

    #ifdef GPTL
        GPTLstart("----mem_split_yx");
    #endif

    mem_split_yx_real(src, src_size[0], src_size[1], src_size[2], work1_r, dims[0], decomp->y1st, decomp->y1en, decomp->y1disp);

    #ifdef GPTL
        GPTLstop("----mem_split_yx");
        GPTLstart("----mpi_alltoall_yx");
    #endif

    MPI_Alltoallv(work1_r, decomp->y1cnts, decomp->y1disp, DECOMP_2D_MPI_REAL,
                  work2_r, decomp->x1cnts, decomp->x1disp, DECOMP_2D_MPI_REAL,
                  decomp_2d_comm_col);

    #ifdef GPTL
        GPTLstop("----mpi_alltoall_yx");
        GPTLstart("----mem_merge_yx");
    #endif

    mem_merge_yx_real(work2_r, dst_size[0], dst_size[1], dst_size[2], dst, dims[0], decomp->x1st, decomp->x1en, decomp->x1disp);

    #ifdef GPTL
        GPTLstop("----mem_merge_yx");
    #endif
}

// Transpose data (real) from Y to Z pencil
void transpose_y_to_z_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp) {

    #ifdef GPTL
        GPTLstart("----mem_split_yz");
    #endif

    mem_split_yz_real(src, src_size[0], src_size[1], src_size[2], work1_r, dims[1], decomp->y2st, decomp->y2en, decomp->y2disp);

    #ifdef GPTL
        GPTLstop("----mem_split_yz");
        GPTLstart("----mpi_alltoall_yz");
    #endif

    MPI_Alltoallv(work1_r, decomp->y2cnts, decomp->y2disp, DECOMP_2D_MPI_REAL,
                  dst,     decomp->z2cnts, decomp->z2disp, DECOMP_2D_MPI_REAL,
                  decomp_2d_comm_row);

    #ifdef GPTL
        GPTLstop("----mpi_alltoall_yz");
        GPTLstart("----mem_merge_yz");
    #endif

    // note the receive buffer is already in natural (i,j,k) order
    // so no merge operation needed

    #ifdef GPTL
        GPTLstop("----mem_merge_yz");
    #endif
}

// Transpose data (real) from Z to Y pencil
void transpose_z_to_y_real(DECOMP_2D_REAL *src, DECOMP_2D_REAL *dst, int *src_size, int *dst_size, decomp_2d_info *decomp) {

    #ifdef GPTL
        GPTLstart("----mem_split_zy");
    #endif

    // note the src array is suitable to be a send buffer
    // so no split operation needed

    #ifdef GPTL
        GPTLstop("----mem_split_zy");
        GPTLstart("----mpi_alltoall_yz");
    #endif

    MPI_Alltoallv(src,     decomp->z2cnts, decomp->z2disp, DECOMP_2D_MPI_REAL,
                  work2_r, decomp->y2cnts, decomp->y2disp, DECOMP_2D_MPI_REAL,
                  decomp_2d_comm_row);

    #ifdef GPTL
        GPTLstop("----mpi_alltoall_yz");
        GPTLstart("----mem_merge_yz");
    #endif

    mem_merge_zy_real(work2_r, dst_size[0], dst_size[1], dst_size[2], dst, dims[1], decomp->y2st, decomp->y2en, decomp->y2disp);

    #ifdef GPTL
        GPTLstop("----mem_merge_yz");
    #endif
}
