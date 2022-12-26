#define UNROLL_F 4
// #define FFT_MKL_CESTOR_CC

#ifdef FFT_MKL_CESTOR_CC
#define DFTI_PLAN dfti_plan_cestorage_cc
#else
#define DFTI_PLAN dfti_plan_cestorage_cr
#endif

// FFTW KIND constants borrowed from the header file fftw3.h
#define FFTW_FORWARD (-1)
#define FFTW_BACKWARD (+1)
enum fftw_r2r_kind_do_not_use_me {
  FFTW_R2HC = 0,
  FFTW_HC2R = 1,
  FFTW_DHT = 2,
  FFTW_REDFT00 = 3,
  FFTW_REDFT01 = 4,
  FFTW_REDFT10 = 5,
  FFTW_REDFT11 = 6,
  FFTW_RODFT00 = 7,
  FFTW_RODFT01 = 8,
  FFTW_RODFT10 = 9,
  FFTW_RODFT11 = 10
};

enum mkl_packed_fmt_udf {
  MKL_PACKED_FMT_PACK = 11,
  MKL_PACKED_FMT_PERM = 12,
  MKL_PACKED_FMT_PERM2 = 13
};

const int mkl_packed_fmt = MKL_PACKED_FMT_PERM;

static void dfti_fwd_posproc_pack(fft_mkl_plan_desc *desc, double *work) {
  memcpy(work + 1, work + 2, desc->n * sizeof(double));
}

static void dfti_fwd_posproc_perm(fft_mkl_plan_desc *desc, double *work) {
  if (desc->n % 2 == 0) {
    work[1] = work[desc->n];
  } else {
    memcpy(work + 1, work + 2, desc->n * sizeof(double));
  }
}

static void dfti_fwd_posproc_perm2(fft_mkl_plan_desc *desc, double *work) {
  work[1] = work[desc->n];
}

static void dfti_bwd_preproc_pack(fft_mkl_plan_desc *desc, double *work) {
  memcpy(work + 2, work + 1, desc->n * sizeof(double));
  work[1] = 0.0;
  if (desc->n % 2 == 0)
    work[desc->n + 1] = 0.0;
}

static void dfti_bwd_preproc_perm(fft_mkl_plan_desc *desc, double *work) {
  if (desc->n % 2 == 0) {
    work[desc->n] = work[1];
    work[desc->n + 1] = 0.0;
  } else {
    memcpy(work + 2, work + 1, desc->n * sizeof(double));
  }
  work[1] = 0.0;
}

static void dfti_bwd_preproc_perm2(fft_mkl_plan_desc *desc, double *work) {
  work[desc->n] = work[1];
  work[desc->n + 1] = 0.0;
  work[1] = 0.0;
}

static void dfti_execute_fi(fft_mkl_plan_desc *desc, double *work) {
  MKL_LONG status = DftiComputeForward(desc->dfti_desc, work);
}

static void dfti_execute_bi(fft_mkl_plan_desc *desc, double *work) {
  MKL_LONG status = DftiComputeBackward(desc->dfti_desc, work);
}

static void dfti_destroy(fft_mkl_plan_desc *desc) {
  MKL_LONG status = DftiFreeDescriptor(&desc->dfti_desc);

  if (desc->ipar != NULL)
    free(desc->ipar);
  if (desc->dpar != NULL)
    free(desc->dpar);
  if (desc->buffer != NULL)
    free(desc->buffer);
}

// static void tt_preproc(fft_mkl_plan_desc *desc, double *work) {}
// static void tt_posproc(fft_mkl_plan_desc *desc, double *work) {}

static void tt_execute_f(fft_mkl_plan_desc *desc, double *work) {
  int ir;
  d_forward_trig_transform(work, &desc->dfti_desc, desc->ipar, desc->dpar, &ir);
}

static void tt_execute_b(fft_mkl_plan_desc *desc, double *work) {
  int ir;
  d_backward_trig_transform(work, &desc->dfti_desc, desc->ipar, desc->dpar,
                            &ir);
}

static void tt_destroy(fft_mkl_plan_desc *desc) {
  int ir;
  free_trig_transform(&desc->dfti_desc, desc->ipar, &ir);

  if (desc->ipar != NULL)
    free(desc->ipar);
  if (desc->dpar != NULL)
    free(desc->dpar);
  if (desc->buffer != NULL)
    free(desc->buffer);
}

static fft_mkl_plan_ptr dfti_plan_cestorage_cc(int fwd_or_bwd, int packed_fmt,
                                               int n, int howmany, int istride,
                                               int idist, int ostride,
                                               int odist) {
  MKL_LONG status, strides[2];

  fft_mkl_plan_ptr plan = (fft_mkl_plan_ptr)malloc(sizeof(fft_mkl_plan));
  plan->desc = (fft_mkl_plan_desc *)malloc(sizeof(fft_mkl_plan_desc));

  plan->desc->packed_fmt = packed_fmt;
  plan->desc->fft_type = 2;
  plan->desc->n = n;
  plan->desc->howmany = howmany;
  plan->desc->istride = istride;
  plan->desc->ostride = ostride;
  plan->desc->idist = idist;
  plan->desc->odist = odist;
  plan->desc->ipar = NULL;
  plan->desc->dpar = NULL;
  plan->desc->buffer = (double *)malloc((n + 2) * UNROLL_F * sizeof(double));

  status = DftiCreateDescriptor(&plan->desc->dfti_desc, DFTI_DOUBLE, DFTI_REAL,
                                1, n);
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_CONJUGATE_EVEN_STORAGE,
                        DFTI_COMPLEX_COMPLEX);
  status =
      DftiSetValue(plan->desc->dfti_desc, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_NUMBER_OF_TRANSFORMS, 1);
  strides[0] = 0;
  strides[1] = 1;
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_INPUT_STRIDES, strides);
  strides[1] = 1;
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_OUTPUT_STRIDES, strides);
  status = DftiCommitDescriptor(plan->desc->dfti_desc);

  if (fwd_or_bwd == FFTW_FORWARD) {
    plan->execute = dfti_execute_fi;
    plan->preproc = NULL;

    switch (packed_fmt) {
    case MKL_PACKED_FMT_PACK:
      plan->posproc = dfti_fwd_posproc_pack;
      break;
    case MKL_PACKED_FMT_PERM:
      plan->posproc = dfti_fwd_posproc_perm;
      break;
    case MKL_PACKED_FMT_PERM2:
      plan->posproc = dfti_fwd_posproc_perm2;
      break;
    }
  } else if (fwd_or_bwd == FFTW_BACKWARD) {
    plan->execute = dfti_execute_bi;
    plan->posproc = NULL;

    switch (packed_fmt) {
    case MKL_PACKED_FMT_PACK:
      plan->preproc = dfti_bwd_preproc_pack;
      break;
    case MKL_PACKED_FMT_PERM:
      plan->preproc = dfti_bwd_preproc_perm;
      break;
    case MKL_PACKED_FMT_PERM2:
      plan->preproc = dfti_bwd_preproc_perm2;
      break;
    }
  }

  plan->destroy = dfti_destroy;

  return plan;
}

static fft_mkl_plan_ptr dfti_plan_cestorage_cr(int fwd_or_bwd, int packed_fmt,
                                               int n, int howmany, int istride,
                                               int idist, int ostride,
                                               int odist) {
  int fmt;
  MKL_LONG status, strides[2]; // NOTE: must be MKL_LONG!!!

  fft_mkl_plan_ptr plan = (fft_mkl_plan_ptr)malloc(sizeof(fft_mkl_plan));
  plan->desc = (fft_mkl_plan_desc *)malloc(sizeof(fft_mkl_plan_desc));

  switch (packed_fmt) {
  case MKL_PACKED_FMT_PACK:
    plan->desc->packed_fmt = packed_fmt;
    fmt = DFTI_PACK_FORMAT;
    break;
  case MKL_PACKED_FMT_PERM:
    plan->desc->packed_fmt = packed_fmt;
    fmt = DFTI_PERM_FORMAT;
    break;
  case MKL_PACKED_FMT_PERM2:
    // MKL does not support PERM2 format
    plan->desc->packed_fmt = MKL_PACKED_FMT_PERM;
    fmt = DFTI_PERM_FORMAT;
    break;
  }
  plan->desc->fft_type = 1;
  plan->desc->n = n;
  plan->desc->howmany = howmany;
  plan->desc->istride = istride;
  plan->desc->ostride = ostride;
  plan->desc->idist = idist;
  plan->desc->odist = odist;
  plan->desc->ipar = NULL;
  plan->desc->dpar = NULL;
  plan->desc->buffer = NULL;

  status = DftiCreateDescriptor(&plan->desc->dfti_desc, DFTI_DOUBLE, DFTI_REAL,
                                1, n);
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_CONJUGATE_EVEN_STORAGE,
                        DFTI_COMPLEX_REAL);
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_PACKED_FORMAT, fmt);
  status =
      DftiSetValue(plan->desc->dfti_desc, DFTI_NUMBER_OF_TRANSFORMS, howmany);
  strides[0] = 0;
  strides[1] = istride;
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_INPUT_STRIDES, strides);
  strides[1] = ostride;
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_OUTPUT_STRIDES, strides);
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_INPUT_DISTANCE, idist);
  status = DftiSetValue(plan->desc->dfti_desc, DFTI_OUTPUT_DISTANCE, odist);
  status = DftiCommitDescriptor(plan->desc->dfti_desc);

  plan->preproc = NULL;
  plan->posproc = NULL;

  if (fwd_or_bwd == FFTW_FORWARD) {
    plan->execute = dfti_execute_fi;
  } else if (fwd_or_bwd == FFTW_BACKWARD) {
    plan->execute = dfti_execute_bi;
  }

  plan->destroy = dfti_destroy;

  return plan;
}

static fft_mkl_plan_ptr tt_plan(int tt_type, int fwd_or_bwd, int n, int howmany,
                                int istride, int idist, int ostride,
                                int odist) {
  fft_mkl_plan_ptr plan = (fft_mkl_plan_ptr)malloc(sizeof(fft_mkl_plan));
  plan->desc = (fft_mkl_plan_desc *)malloc(sizeof(fft_mkl_plan_desc));

  plan->desc->packed_fmt = -1; // unused
  plan->desc->fft_type = 3;
  plan->desc->n = n;
  plan->desc->howmany = howmany;
  plan->desc->istride = istride;
  plan->desc->ostride = ostride;
  plan->desc->idist = idist;
  plan->desc->odist = odist;

  plan->desc->ipar = (MKL_INT *)malloc(128 * sizeof(MKL_INT));
  plan->desc->dpar = (double *)malloc((5 * n / 2 + 2) * sizeof(double));

  int ir;
  d_init_trig_transform(&n, &tt_type, plan->desc->ipar, plan->desc->dpar, &ir);
  plan->desc->ipar[10] = 1; // compatibility with FFTW, unnormalized
  double *f = (double *)malloc(n * sizeof(double));
  d_commit_trig_transform(f, &plan->desc->dfti_desc, plan->desc->ipar,
                          plan->desc->dpar, &ir);
  free(f);

  if (istride > 1) {
    plan->desc->buffer = (double *)malloc(UNROLL_F * n);
  } else {
    plan->desc->buffer = NULL;
  }

  plan->preproc = NULL; // not yet implemented for DFTI_COMPLEX_COMPLEX FFT
  plan->posproc = NULL; // not yet implemented for DFTI_COMPLEX_COMPLEX FFT

  if (fwd_or_bwd == FFTW_FORWARD) {
    plan->execute = tt_execute_f;
  } else if (fwd_or_bwd == FFTW_BACKWARD) {
    plan->execute = tt_execute_b;
  }

  plan->destroy = tt_destroy;

  return plan;
}

// Workarounds, not equivalent to FFTW, but useful for solving PDE by
// eigenfunction expansion method. Note that the vector containing corresponding
// eigen values should be reorganized to match with the MKL packed format.
static fft_mkl_plan_ptr fft_mkl_plan_many_r2r_1d(const int packed_fmt,
                                                 const int fftw_kind, int n,
                                                 int howmany, int istride,
                                                 int idist, int ostride,
                                                 int odist) {
  switch (fftw_kind) {

  case FFTW_DHT:
    return DFTI_PLAN(FFTW_FORWARD, packed_fmt, n, howmany, istride, idist,
                     ostride, odist);
    break;

  case FFTW_R2HC:
    return DFTI_PLAN(FFTW_FORWARD, packed_fmt, n, howmany, istride, idist,
                     ostride, odist);
    break;

  case FFTW_HC2R:
    return DFTI_PLAN(FFTW_BACKWARD, packed_fmt, n, howmany, istride, idist,
                     ostride, odist);
    break;

  case FFTW_REDFT10:
    return tt_plan(MKL_STAGGERED_COSINE_TRANSFORM, FFTW_BACKWARD, n, howmany,
                   istride, idist, ostride, odist);
    break;

  case FFTW_REDFT01:
    return tt_plan(MKL_STAGGERED_COSINE_TRANSFORM, FFTW_FORWARD, n, howmany,
                   istride, idist, ostride, odist);
    break;

  case FFTW_REDFT11:
    return tt_plan(MKL_STAGGERED2_COSINE_TRANSFORM, FFTW_BACKWARD, n, howmany,
                   istride, idist, ostride, odist);
    break;

  case FFTW_RODFT10:
    return tt_plan(MKL_STAGGERED_SINE_TRANSFORM, FFTW_BACKWARD, n, howmany,
                   istride, idist, ostride, odist);
    break;

  case FFTW_RODFT01:
    return tt_plan(MKL_STAGGERED_SINE_TRANSFORM, FFTW_FORWARD, n, howmany,
                   istride, idist, ostride, odist);
    break;

  case FFTW_RODFT11:
    return tt_plan(MKL_STAGGERED2_SINE_TRANSFORM, FFTW_BACKWARD, n, howmany,
                   istride, idist, ostride, odist);
    break;
  }
}

static void gather(const int unroll_factor, const int n, const int stride,
                   const int dist, const double *src, const int lead_dim,
                   double *dst) {
  if (stride == 1) {
    for (int j = 0; j < unroll_factor; j++) {
      memcpy(dst + lead_dim * j, src + dist * j, n * sizeof(double));
    }
  } else {
    int pos = 0;
    for (int i = 0; i < n; i++) {
// #pragma GCC unroll 4
#pragma unroll(UNROLL_F)
      for (int j = 0; j < unroll_factor; j++) {
        dst[i + lead_dim * j] = src[pos + dist * j];
      }
      pos += stride;
    }
  }
}

static void scatter(const int unroll_factor, const int lead_dim,
                    const double *src, const int n, const int stride,
                    const int dist, double *dst) {
  if (stride == 1) {
    for (int j = 0; j < unroll_factor; j++) {
      memcpy(dst + dist * j, src + lead_dim * j, n * sizeof(double));
    }
  } else {
    int pos = 0;
    for (int i = 0; i < n; i++) {
// #pragma GCC unroll 4
#pragma unroll(UNROLL_F)
      for (int j = 0; j < unroll_factor; j++) {
        dst[pos + dist * j] = src[i + lead_dim * j];
      }
      pos += stride;
    }
  }
}

static void fft_mkl_execute_r2r(fft_mkl_plan_ptr plan, double *work) {

  // execute in-place FFT
  switch (plan->desc->fft_type) {
  case 1:
    // dfti in DFTI_COMPLEX_REAL storage
    plan->execute(plan->desc, work);
    break;

  case 2: {
    // dfti in DFTI_COMPLEX_COMPLEX storage
    int buf_lead_dim = plan->desc->n + 2;
    int irest = plan->desc->howmany % UNROLL_F;
    int pos;
    int i = 0;
    for (; i < irest; i++) {
      pos = i * plan->desc->idist;
      gather(1, plan->desc->n, plan->desc->istride, plan->desc->idist,
             work + pos, buf_lead_dim, plan->desc->buffer);
      if (plan->preproc != NULL)
        plan->preproc(plan->desc, plan->desc->buffer);
      plan->execute(plan->desc, plan->desc->buffer);
      if (plan->posproc != NULL)
        plan->posproc(plan->desc, plan->desc->buffer);
      scatter(1, buf_lead_dim, plan->desc->buffer, plan->desc->n,
              plan->desc->istride, plan->desc->idist, work + pos);
    }
    for (; i < plan->desc->howmany; i += UNROLL_F) {
      pos = i * plan->desc->idist;
      gather(UNROLL_F, plan->desc->n, plan->desc->istride, plan->desc->idist,
             work + pos, buf_lead_dim, plan->desc->buffer);
// #pragma GCC unroll 4
#pragma unroll(UNROLL_F)
      for (int j = 0; j < UNROLL_F; j++) {
        int iloc = buf_lead_dim * j;
        if (plan->preproc != NULL)
          plan->preproc(plan->desc, plan->desc->buffer + iloc);
        plan->execute(plan->desc, plan->desc->buffer + iloc);
        if (plan->posproc != NULL)
          plan->posproc(plan->desc, plan->desc->buffer + iloc);
      }
      scatter(UNROLL_F, buf_lead_dim, plan->desc->buffer, plan->desc->n,
              plan->desc->istride, plan->desc->idist, work + pos);
    }
    break;
  }
  case 3:
    // tt
    if (plan->desc->istride > 1) {
      int buf_lead_dim = plan->desc->n + 2;
      int irest = plan->desc->howmany % UNROLL_F;
      int pos;
      int i = 0;
      for (; i < irest; i++) {
        pos = i * plan->desc->idist;
        gather(1, plan->desc->n, plan->desc->istride, plan->desc->idist,
               work + pos, buf_lead_dim, plan->desc->buffer);
        if (plan->preproc != NULL)
          plan->preproc(plan->desc, plan->desc->buffer);
        plan->execute(plan->desc, plan->desc->buffer);
        if (plan->posproc != NULL)
          plan->posproc(plan->desc, plan->desc->buffer);
        scatter(1, buf_lead_dim, plan->desc->buffer, plan->desc->n,
                plan->desc->istride, plan->desc->idist, work + pos);
      }
      for (; i < plan->desc->howmany; i += UNROLL_F) {
        pos = i * plan->desc->idist;
        gather(UNROLL_F, plan->desc->n, plan->desc->istride, plan->desc->idist,
               work + pos, buf_lead_dim, plan->desc->buffer);
// #pragma GCC unroll 4
#pragma unroll(UNROLL_F)
        for (int j = 0; j < UNROLL_F; j++) {
          int iloc = buf_lead_dim * j;
          if (plan->preproc != NULL)
            plan->preproc(plan->desc, plan->desc->buffer + iloc);
          plan->execute(plan->desc, plan->desc->buffer + iloc);
          if (plan->posproc != NULL)
            plan->posproc(plan->desc, plan->desc->buffer + iloc);
        }
        scatter(UNROLL_F, buf_lead_dim, plan->desc->buffer, plan->desc->n,
                plan->desc->istride, plan->desc->idist, work + pos);
      }
    } else {
      for (int i = 0; i < plan->desc->howmany; i++) {
        int pos = i * plan->desc->idist;
        if (plan->preproc != NULL)
          plan->preproc(plan->desc, work + pos);
        plan->execute(plan->desc, work + pos);
        if (plan->posproc != NULL)
          plan->posproc(plan->desc, work + pos);
      }
    }
    break;
  }
}

static void fft_mkl_destroy_plan(fft_mkl_plan_ptr plan) {
  plan->destroy(plan->desc);
  if (plan->desc != NULL)
    free(plan->desc);
  if (plan != NULL)
    free(plan);
}

static void get_fftw_kind(char bctype[], int *kind_fwd, int *kind_bwd) {
  if (strcmp(bctype, "PP") == 0) {
    *kind_fwd = FFTW_R2HC;
    *kind_bwd = FFTW_HC2R;
  } else if (strcmp(bctype, "NN") == 0) {
    *kind_fwd = FFTW_REDFT10;
    *kind_bwd = FFTW_REDFT01;
  } else if (strcmp(bctype, "DD") == 0) {
    *kind_fwd = FFTW_RODFT10;
    *kind_bwd = FFTW_RODFT01;
  } else if (strcmp(bctype, "ND") == 0) {
    *kind_fwd = FFTW_REDFT11;
    *kind_bwd = FFTW_REDFT11;
  } else if (strcmp(bctype, "DN") == 0) {
    *kind_fwd = FFTW_RODFT11;
    *kind_bwd = FFTW_RODFT11;
  }
}

static void get_norm_factor(char bctype[], double *normfactor) {
  if (strcmp(bctype, "PP") == 0) {
    *normfactor = 1.0;
  } else if (strcmp(bctype, "NN") == 0) {
    *normfactor = 2.0;
  } else if (strcmp(bctype, "DD") == 0) {
    *normfactor = 2.0;
  } else if (strcmp(bctype, "ND") == 0) {
    *normfactor = 2.0;
  } else if (strcmp(bctype, "DN") == 0) {
    *normfactor = 2.0;
  }
}

void init_fft(int xsz[3], int ysz[3], char bctype_x[], char bctype_y[],
              fft_mkl_plan_ptr fft_plan[2][2], double *fft_normfactor) {
  fft_mkl_plan_ptr plan_fwd_xpen, plan_bwd_xpen;
  fft_mkl_plan_ptr plan_fwd_ypen, plan_bwd_ypen;

  int kind_fwd, kind_bwd;

  // Init single in-place FFT along x direction (x-pencil)
  double normfactor_x;
  get_fftw_kind(bctype_x, &kind_fwd, &kind_bwd);
  get_norm_factor(bctype_x, &normfactor_x);
  plan_fwd_xpen = fft_mkl_plan_many_r2r_1d(mkl_packed_fmt, kind_fwd, xsz[0], xsz[1],
                                           1, xsz[0], 1, xsz[0]);
  plan_bwd_xpen = fft_mkl_plan_many_r2r_1d(mkl_packed_fmt, kind_bwd, xsz[0], xsz[1],
                                           1, xsz[0], 1, xsz[0]);

  // Init single in-place FFT along y direction (y-pencil)
  double normfactor_y;
  get_fftw_kind(bctype_y, &kind_fwd, &kind_bwd);
  get_norm_factor(bctype_y, &normfactor_y);
  plan_fwd_ypen = fft_mkl_plan_many_r2r_1d(mkl_packed_fmt, kind_fwd, ysz[1], ysz[0],
                                           ysz[0], 1, ysz[0], 1);
  plan_bwd_ypen = fft_mkl_plan_many_r2r_1d(mkl_packed_fmt, kind_bwd, ysz[1], ysz[0],
                                           ysz[0], 1, ysz[0], 1);
  // plan_fwd_ypen = fft_mkl_plan_many_r2r_1d(mkl_packed_fmt, kind_fwd, ysz[1],
  // 1, 1, ysz[0], 1, ysz[0]);
  // plan_bwd_ypen = fft_mkl_plan_many_r2r_1d(mkl_packed_fmt, kind_fwd, ysz[1],
  // 1, 1, ysz[0], 1, ysz[0]);

  fft_plan[0][0] = plan_fwd_xpen;
  fft_plan[0][1] = plan_bwd_xpen;
  fft_plan[1][0] = plan_fwd_ypen;
  fft_plan[1][1] = plan_bwd_ypen;
  *fft_normfactor = 1.0 / (normfactor_x * xsz[0] * normfactor_y * ysz[1]);
}

void execute_fft(fft_mkl_plan_ptr plan, double *work) {
  // Execute in-place FFT
  fft_mkl_execute_r2r(plan, work);
}

void free_fft(fft_mkl_plan_ptr plan[2][2]) {
  for (int idir = 0; idir < 2; idir++) {
    for (int i = 0; i < 2; i++) {
      fft_mkl_destroy_plan(plan[idir][i]);
    }
  }
}

void get_eigen_values(int ist, int isz, int isz_global, char bctype[],
                      double *lambda) {
  const double pi = acos(-1.0);
  int ien = ist + isz - 1;

  if (strcmp(bctype, "PP") == 0) {
    int n = isz_global;
    double *lambda_glb = (double *)malloc(n * sizeof(double));
    double *lambda_aux = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
      lambda_aux[i] = 2.0 * (cos(2.0 * pi * (i - 0.0) / n) - 1.0);
    }

    switch (mkl_packed_fmt) {
    case MKL_PACKED_FMT_PACK:
      lambda_glb[0] = lambda_aux[0];
      for (int i = 1; i < n; i++) {
        lambda_glb[i] = lambda_aux[(i + 1) / 2];
      }
      break;

    case MKL_PACKED_FMT_PERM:
      if (n % 2 == 0) {
        lambda_glb[0] = lambda_aux[0];
        lambda_glb[1] = lambda_aux[n / 2];
        for (int i = 2; i < n; i++) {
          lambda_glb[i] = lambda_aux[i / 2];
        }
      } else {
        lambda_glb[0] = lambda_aux[0];
        for (int i = 1; i < n; i++) {
          lambda_glb[i] = lambda_aux[(i + 1) / 2];
        }
      }
      break;

    case MKL_PACKED_FMT_PERM2:
      lambda_glb[0] = lambda_aux[0];
      lambda_glb[1] = lambda_aux[n / 2];
      for (int i = 2; i < n; i++) {
        lambda_glb[i] = lambda_aux[i / 2];
      }
      break;
    }

    for (int i = ist; i <= ien; i++) {
      lambda[i - ist] = lambda_glb[i - 1];
    }

    free(lambda_glb);
    free(lambda_aux);

  } else if (strcmp(bctype, "NN") == 0) {
    for (int i = ist; i <= ien; i++) {
      lambda[i - ist] = 2.0 * (cos(pi * (i - 1.0) / isz_global) - 1.0);
    }
  } else if (strcmp(bctype, "DD") == 0) {
    for (int i = ist; i <= ien; i++) {
      lambda[i - ist] = 2.0 * (cos(pi * (i - 0.0) / isz_global) - 1.0);
    }
  } else if (strcmp(bctype, "ND") == 0 || strcmp(bctype, "DN") == 0) {
    for (int i = ist; i <= ien; i++) {
      lambda[i - ist] =
          2.0 * (cos(pi * (2 * i - 1.0) / (2.0 * isz_global)) - 1.0);
    }
  }
}
