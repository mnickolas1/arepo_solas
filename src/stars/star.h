#ifndef STAR_H
#define STAR_H

#include "../main/allvars.h"

#ifdef STAR_PARTICLES
#include "../stars/star_particle.h"
#endif

#ifdef STAR_RADIATION_ACTIVE
#include "../stars/star_radiation.h"
#endif

#ifdef STAR_FEEDBACK_ACTIVE
#include "../stars/star_tables.h"
#endif


#define ALLOC_STAR_ROOM 64
extern int NumStars;

#ifdef STAR_FEEDBACK_ACTIVE
extern struct TimeBinData TimeBinsStar;

typedef struct Star_Interpolate
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
  WavebandData Flux[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Star_Interpolate;

typedef struct Star_Feedback
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
  WavebandData Radiated[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Star_Feedback;

typedef struct Mechanical_Feedback
{
#ifdef WINDS 
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE
  WavebandData Radiated[WAVEBANDS];
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif
} Mechanical_Feedback;
#endif

typedef struct Star_Particle_Data
{
  MyIDType PID;

#ifdef METALS
  MyDouble Metallicity;
#endif

#ifdef STAR_PARTICLES
  MyDouble MassOfStar;
#endif

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
  int NumOfStarsInBins[NBINS];
#endif

#ifdef INDIVIDUAL_STAR_BY_STAR_FORMATION
  MyDouble MassToDrain;
#endif

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeSN_yr;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
  int Active;
  int WithFeedback;
  MyDouble Hsml;
  MyDouble NgbsMass;
  MyDouble NgbsVolume;
#ifdef METALS
  MyDouble NgbsMetallicity;
#endif
  int NgbsMinBin;
  int DensityFlag;
  signed char TimeBinStar;
  MyDouble PhysicalAge_yr;
  Mechanical_Feedback MechanicalFeedback;
#endif
} Star_Particle_Data;

extern Star_Particle_Data *SP;

#define SPP(i) SP[P[i].SID]
#define PPS(i) P[SP[i].PID]

#endif