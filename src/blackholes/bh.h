#ifndef BH_H
#define BH_H

#include "../utils/dtypes.h"

#define ALLOC_BH_ROOM 4
extern int NumBhs;

#if defined(BH_ACCRETION_ACTIVE) || defined(BH_FEEDBACK_ACTIVE)
#include "../time_integration/timestep.h"
extern struct TimeBinData TimeBinsBh;
#endif

extern struct bh_particle_data
{
  MyIDType PID;
  MyDouble Hsml;
  MyDouble Density;
  MyDouble NgbMass;
  MyDouble NgbMassFeed;
#ifdef BONDI_ACCRETION
  MyDouble VelocityGas[3];
  MyDouble VelocityGasCircular[3];
  MyDouble InternalEnergyGas;
  MyDouble AccretionRate;
  MyDouble MassToDrain;
  MyDouble AngularMomentum[3];
#endif
#ifdef INFALL_ACCRETION
  MyDouble Accretion;
#endif
  integertime NgbMinStep;
  int DensityFlag;
  signed char TimeBinBh;
}  *BhP;

#define BPP(i) BhP[P[i].BhID]
#define PPB(i) P[BhP[i].PID]

#endif