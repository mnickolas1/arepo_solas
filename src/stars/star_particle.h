#ifndef STAR_PARTICLE_H
#define STAR_PARTICLE_H

#define NBINS 114
#define MMIN 0.1
#define MMAX 120

#define N_CDF_BINS 10000

extern double cdf_masses[N_CDF_BINS + 1];   
extern double cdf_values[N_CDF_BINS + 1];   

#if STAR_PARTICLES == 1
extern double StarMassBins[NBINS + 1];
extern double StarMeanMassInBins[NBINS];
#endif

#endif