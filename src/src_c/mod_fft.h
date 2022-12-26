// #define USE_MKL
#ifdef USE_MKL

#include "mkl_dfti.h"
#include "mkl_trig_transforms.h"

struct FFT_MKL_PLAN_DESC {
  DFTI_DESCRIPTOR_HANDLE dfti_desc;
  int fft_type;   // 1 - dfti in DFTI_COMPLEX_REAL storage; 2 - dfti in
                  // DFTI_COMPLEX_COMPLEX storage; 3 - tt
  int packed_fmt; // 11 - DFTI_PACK_FORMAT; 12 - DFTI_PERM_FORMAT; 13 - PERM2
                  // format (new); -1 - unused
  int n, howmany;
  int istride, ostride;
  int idist, odist;
  MKL_INT *ipar;
  double *dpar;
  double *buffer;
};
typedef struct FFT_MKL_PLAN_DESC fft_mkl_plan_desc;

typedef void (*prepost_funptr)(fft_mkl_plan_desc *, double *);
typedef void (*execute_funptr)(fft_mkl_plan_desc *, double *);
typedef void (*destroy_funptr)(fft_mkl_plan_desc *);
typedef struct FFT_MKL_PLAN fft_mkl_plan;
typedef fft_mkl_plan *fft_mkl_plan_ptr;
struct FFT_MKL_PLAN {
  fft_mkl_plan_desc *desc;
  prepost_funptr preproc;
  execute_funptr execute;
  prepost_funptr posproc;
  destroy_funptr destroy;
};

extern void init_fft(int xsz[3], int ysz[3], char bctype_x[], char bctype_y[],
                     fft_mkl_plan_ptr fft_plan[2][2], double *fft_normfactor);
extern void execute_fft(fft_mkl_plan_ptr plan, double *work);
extern void free_fft(fft_mkl_plan_ptr plan[2][2]);
extern void get_eigen_values(int ist, int isz, int isz_global, char bctype[],
                             double *lambda);

#else

#include <fftw3.h>

extern void init_fft(int xsz[3], int ysz[3], char bctype_x[], char bctype_y[],
                     double *work_xpen, double *work_ypen,
                     fftw_plan fft_plan[2][2], double *fft_normfactor);
extern void execute_fft(fftw_plan plan, double *work);
// extern void execute_fft(fftw_plan plan, double *in, double *out);
extern void free_fft(fftw_plan plan[2][2]);
extern void get_eigen_values(int ist, int isz, int isz_global, char bctype[],
                             double *lambda);

#endif
