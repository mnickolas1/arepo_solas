#ifndef STAR_TABLES_H
#define STAR_TABLES_H

#include "../main/allvars.h"


extern int Z_COUNT;
extern int M_COUNT;

extern double *Z_VALUES;
extern double *M_VALUES;

extern double *logZ_VALUES;
extern double *logM_VALUES;

extern int **N;

extern double ***Age;
extern double ***FractionalAge;
extern double ***logRadius;
extern double ***logTemperature;

#ifdef WINDS
extern double ***logMassLossRate;
#if GRACKLE_CHEMISTRY >= 1
extern double ***WindX;
extern double ***WindY;
#endif
#ifdef METALS
extern double ***WindZ;
#endif
extern double ***logWindVelocity;
#endif

#ifdef STAR_RADIATION_ACTIVE
extern WavebandData ***logFlux[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
extern double **SN_MassLoss; 
#if GRACKLE_CHEMISTRY >= 1
extern double **SN_X; 
extern double **SN_Y; 
#endif
#ifdef METALS
extern double **SN_Z; 
#endif 
#endif

#ifdef AGB 
extern double **AGB_MassLoss; 
#ifdef METALS
extern double **AGB_MetalsLoss; 
#endif 
#endif

#endif