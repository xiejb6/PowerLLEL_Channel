static void get_fftw_kind(char bctype[], fftw_r2r_kind *kind_fwd,
                          fftw_r2r_kind *kind_bwd) {
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
              double *work_xpen, double *work_ypen, fftw_plan fft_plan[2][2],
              double *fft_normfactor) {
  fftw_plan plan_fwd_xpen, plan_bwd_xpen;
  fftw_plan plan_fwd_ypen, plan_bwd_ypen;

  const int rank = 1;
  const int howmany_rank = 1;
  fftw_iodim dims_t[rank];
  fftw_iodim howmany_dims_t[howmany_rank];
  fftw_r2r_kind kind_fwd[rank], kind_bwd[rank];

  // Init single in-place FFT along x direction (x-pencil)
  dims_t[0].n = xsz[0];
  dims_t[0].is = 1;
  dims_t[0].os = 1;
  howmany_dims_t[0].n = xsz[1];
  howmany_dims_t[0].is = xsz[0]; // unused
  howmany_dims_t[0].os = xsz[0]; // unused

  double normfactor_x;
  get_fftw_kind(bctype_x, &kind_fwd[0], &kind_bwd[0]);
  get_norm_factor(bctype_x, &normfactor_x);
  plan_fwd_xpen =
      fftw_plan_guru_r2r(rank, dims_t, howmany_rank, howmany_dims_t, work_xpen,
                         work_xpen, kind_fwd, FFTW_MEASURE);
  plan_bwd_xpen =
      fftw_plan_guru_r2r(rank, dims_t, howmany_rank, howmany_dims_t, work_xpen,
                         work_xpen, kind_bwd, FFTW_MEASURE);

  // Init single in-place FFT along y direction (y-pencil)
  dims_t[0].n = ysz[1];
  // dims_t[0].is = 1;
  // dims_t[0].os = 1;
  dims_t[0].is = ysz[0];
  dims_t[0].os = ysz[0];
  howmany_dims_t[0].n = ysz[0];
  howmany_dims_t[0].is = 1; // unused
  howmany_dims_t[0].os = 1; // unused

  double normfactor_y;
  get_fftw_kind(bctype_y, &kind_fwd[0], &kind_bwd[0]);
  get_norm_factor(bctype_y, &normfactor_y);
  plan_fwd_ypen =
      fftw_plan_guru_r2r(rank, dims_t, howmany_rank, howmany_dims_t, work_ypen,
                         work_ypen, kind_fwd, FFTW_MEASURE);
  plan_bwd_ypen =
      fftw_plan_guru_r2r(rank, dims_t, howmany_rank, howmany_dims_t, work_ypen,
                         work_ypen, kind_bwd, FFTW_MEASURE);

  fft_plan[0][0] = plan_fwd_xpen;
  fft_plan[0][1] = plan_bwd_xpen;
  fft_plan[1][0] = plan_fwd_ypen;
  fft_plan[1][1] = plan_bwd_ypen;
  *fft_normfactor = 1.0 / (normfactor_x * xsz[0] * normfactor_y * ysz[1]);
}

void execute_fft(fftw_plan plan, double *work) {
  // Execute in-place FFT
  fftw_execute_r2r(plan, work, work);
}

// void execute_fft(fftw_plan plan, double *in, double *out) {
//   // Execute out-of-place FFT
//   fftw_execute_r2r(plan, in, out);
// }

void free_fft(fftw_plan plan[2][2]) {
  for (int idir = 0; idir < 2; idir++) {
    for (int i = 0; i < 2; i++) {
      fftw_destroy_plan(plan[idir][i]);
    }
  }
}

void get_eigen_values(int ist, int isz, int isz_global, char bctype[],
                      double *lambda) {
  const double pi = acos(-1.0);
  int ien = ist + isz - 1;

  if (strcmp(bctype, "PP") == 0) {
    for (int i = ist; i <= ien; i++) {
      lambda[i - ist] = 2.0 * (cos(2.0 * pi * (i - 1.0) / isz_global) - 1.0);
    }
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
