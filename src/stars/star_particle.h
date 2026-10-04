#ifndef STAR_PARTICLE_H
#define STAR_PARTICLE_H

#include <limits.h>

#define NBINS 265

#ifndef STAR_BINS_BITS 
#define STAR_BINS_BITS 32
#endif

#if STAR_BINS_BITS == 8
typedef unsigned char MyStarBins;
#define BIN_COUNTS_MAX UCHAR_MAX
#elif STAR_BINS_BITS == 16
typedef unsigned short MyStarBins;
#define BIN_COUNTS_MAX USHRT_MAX
#elif STAR_BINS_BITS == 32
typedef unsigned int MyStarBins;
#define BIN_COUNTS_MAX UINT_MAX
#else
#error "STAR_BINS_BITS must be 8, 16 or 32"
#endif

#define MMIN 0.10
#define MMAX 120.0

#define N_CDF_BINS 10000

extern double cdf_masses[N_CDF_BINS + 1];
extern double cdf_values[N_CDF_BINS + 1];

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
extern double StarMassBins[NBINS + 1];
extern double StarMeanMassInBins[NBINS];
#endif

#if defined(STAR_PARTICLES) && STAR_PARTICLES == 0

#include <gsl/gsl_rng.h>

extern gsl_rng *rng;

extern double norm;
extern double bin_imf[NBINS];
#endif

#endif