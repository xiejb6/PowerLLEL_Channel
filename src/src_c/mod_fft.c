#include "mod_fft.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef USE_MKL
#include "mod_fft_mkl.inc"
#else
#include "mod_fft_fftw.inc"
#endif
