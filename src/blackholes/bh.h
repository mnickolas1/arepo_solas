#ifndef BH_H
#define BH_H

#define ALLOC_BH_ROOM 4
extern int NumBhs;

#include "../main/allvars.h"


#if defined(BH_ACCRETION_ACTIVE) || defined(BH_FEEDBACK_ACTIVE)
extern struct TimeBinData TimeBinsBh;
#endif

extern struct bh_particle_data
{
  MyIDType PID;

#if defined(BH_ACCRETION_ACTIVE) || defined(BH_FEEDBACK_ACTIVE)
  MyDouble Hsml;
  MyDouble NgbMass;
  MyDouble NgbVolume;
  MyDouble AngularMomentum[3];
  int NgbMaxBin;
  int DensityFlag;
  signed char TimeBinBh;
  
  MyDouble NgbMassFeed;
  MyDouble NgbVolumeFeed;
#endif

#ifdef BONDI_ACCRETION
  MyDouble VelocityGas[3];
  MyDouble VelocityGasCircular[3];
  MyDouble Density;
  MyDouble InternalEnergyGas;
  MyDouble AccretionRate;
  MyDouble MassToDrain;
#endif

#ifdef TORQUE_ACCRETION
  MyDouble TorqueMgas;
  MyDouble TorqueMstar;
  MyDouble TorqueMgasDisk;
  MyDouble TorqueMstarDisk;
  MyDouble TorqueR0;
  MyDouble TorqueFd;
  MyDouble VelocityGasCircular[3];
  MyDouble InternalEnergyGas;
  MyDouble AccretionRate;
  MyDouble MassToDrain;
#endif

#ifdef ADP_ACCRETION
  MyDouble ADP_CapturedMass;   /* mass captured since last update (code mass) */
  MyDouble ADP_ReservoirMass;  /* reservoir mass waiting to enter disc (code mass) */
  MyDouble ADP_DiscMass;       /* disc mass available to accrete (code mass) */
  MyDouble VelocityGas[3];
  MyDouble VelocityGasCircular[3];
  MyDouble AccretionRate;
  MyDouble MassToDrain; 
#endif

//#ifdef INFALL_ACCRETION
//  MyDouble Accretion;
//#endif
}  *BhP;

#define BPP(i) BhP[P[i].BhID]
#define PPB(i) P[BhP[i].PID]

#endif