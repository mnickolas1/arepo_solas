#ifndef STAR_TABLES_H
#define STAR_TABLES_H

#include <mpi.h>

extern int Z_COUNT;
extern int M_COUNT;

extern double *Z_VALUES;
extern double *M_VALUES;

extern int **N;

extern double ***Age;
extern double ***Radius;
extern double ***Temperature;

#ifdef WINDS
extern double ***MassLossRate;
#ifdef METALS
extern double ***MetalsLossRate;
#endif
extern double ***WindVelocity;
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
extern double ***RAD_IonizingRate;
extern double ***RAD_IonizingLuminosity;
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
extern double ***RAD_UVLymanWernerLuminosity;
extern double ***RAD_UltravioletLuminosity;
#endif
#if defined(RADIATION_PRESSURE)
extern double ***RAD_OpticalLuminosity;
extern double ***RAD_InfraredLuminosity;
#endif

#ifdef SUPERNOVAE
extern double **SN_MassLoss; 
#ifdef METALS
extern double **SN_MetalsLoss; 
#endif 
#endif

#ifdef AGB 
extern double **AGB_MassLoss; 
#ifdef METALS
extern double **AGB_MetalsLoss; 
#endif 
#endif

#endif