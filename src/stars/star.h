#ifndef STAR_H
#define STAR_H

#define ALLOC_STAR_ROOM 64
extern int NumStars;

#include "../main/allvars.h"

#ifdef STAR_PARTICLES
#include "../stars/star_particle.h"
#endif

#ifdef STAR_FEEDBACK_ACTIVE
#include "../stars/star_tables.h"
#endif

#ifdef STAR_RADIATION_ACTIVE
#include "../stars/star_radiation.h"
#endif


#ifdef STAR_FEEDBACK_ACTIVE
extern struct TimeBinData TimeBinsStar;
#endif

extern struct star_particle_data
{
  MyIDType PID;

#ifdef METALS
  MyDouble Metallicity;
#endif

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
  int NumOfStarsInBins[NBINS];
#endif

#if STAR_PARTICLES == 2
  MyDouble MassOfStar;
#endif

#ifdef INDIVIDUAL_STAR_BY_STAR_FORMATION
  MyDouble MassToDrain;
#endif

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeSN_yr;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
  int Active;
  MyDouble Hsml;
  MyDouble NgbsMass;
  MyDouble NgbsVolume;
  int NgbsMaxBin;
  int DensityFlag;
  signed char TimeBinStar;
  MyDouble PhysicalAge_yr;
#endif

#ifdef WINDS 
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE
  struct WavebandData Radiated[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
}  *SP;

#define SPP(i) SP[P[i].SID]
#define PPS(i) P[SP[i].PID]

#ifdef STAR_FEEDBACK_ACTIVE
struct star_interpolate
{
  MyDouble Radius;
  MyDouble Temperature;

#ifdef WINDS
  MyDouble MassLossRate;
#ifdef METALS
  MyDouble MetalsLossRate;
#endif
  MyDouble WindVelocity;
#endif

#ifdef STAR_RADIATION_ACTIVE
  struct WavebandData Flux[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
};

struct star_feedback
{
#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  double TimeSN;
#endif
  
  int Stage; // 0:preSN, 1:SN, 2:postSN

#ifdef WINDS
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE
  struct WavebandData Radiated[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
};
#endif

#endif