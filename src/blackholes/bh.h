#ifndef BH_H
#define BH_H

#define ALLOC_BH_ROOM 4
extern int NumBhs;

#include "../main/allvars.h"

extern FILE *FdBlackHoles; 

#if defined(BH_ACCRETION_ACTIVE) || defined(BH_FEEDBACK_ACTIVE)
extern struct TimeBinData TimeBinsBh;
#endif

extern struct bh_particle_data
{
  MyIDType PID;

#if defined(BH_ACCRETION_ACTIVE) || defined(BH_FEEDBACK_ACTIVE)
  MyDouble Hsml;
  MyDouble NgbsMass;
  MyDouble NgbsVolume;
  MyDouble AngularMomentum[3];
  int NgbsMaxBin;
  int DensityFlag;
  signed char TimeBinBh;
#endif

#ifdef BH_ACCRETION_ACTIVE
  MyDouble Accretion;
  MyDouble MassToDrain;
#endif

#ifdef BONDI_ACCRETION
  MyDouble GasVelocity[3];
  MyDouble GasCircularVelocity[3];
  MyDouble GasDensity;
  MyDouble GasInternalEnergy;
#endif

#ifdef TORQUE_ACCRETION
  MyDouble TorqueMgas;
  MyDouble TorqueMstar;
  MyDouble TorqueMgasDisk;
  MyDouble TorqueMstarDisk;
  MyDouble TorqueR0;
  MyDouble TorqueFd;
  MyDouble GasAngularMomentum[3];
  MyDouble GasCircularVelocity[3];
#endif

#ifdef ADP_ACCRETION
  MyDouble ADP_CapturedMass;   /* mass captured since last update (code mass) */
  MyDouble ADP_ReservoirMass;  /* reservoir mass waiting to enter disc (code mass) */
  MyDouble ADP_DiscMass;       /* disc mass available to accrete (code mass) */
  MyDouble GasCircularVelocity[3];
#endif

//#ifdef INFALL_ACCRETION
//  MyDouble Accretion;
//#endif
}  *BhP;

#define BPP(i) BhP[P[i].BhID]
#define PPB(i) P[BhP[i].PID]

#endif