#ifndef STAR_H
#define STAR_H

#include "../main/allvars.h"

#define ALLOC_STAR_ROOM 16
extern int NumStars;

#ifdef STAR_FEEDBACK
#define WINDS
#define RADIATION
#define SUPERNOVAE
#endif

#ifdef RADIATION
#define PHOTOIONIZATION
#define PHOTOELECTRIC_HEATING
#define RADIATION_PRESSURE
#endif

#if defined(WINDS) \
|| defined(PHOTOIONIZATION) \
|| defined(PHOTOELECTRIC_HEATING) \
|| defined(RADIATION_PRESSURE) \
|| defined(SUPERNOVAE)
#define STAR_FEEDBACK_ACTIVE
#endif

#ifdef STAR_FEEDBACK_ACTIVE
#ifndef STARS
#error "We cannot run star feeedback simulations without stars!"
#endif
#endif

#if defined(PHOTOIONIZATION) || defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
#define STAR_RADIATION_ACTIVE
#endif

#ifdef STAR_FEEDBACK_ACTIVE
#include "../time_integration/timestep.h"
extern struct TimeBinData TimeBinsStar;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
#ifndef STAR_BY_STAR
#include "../stars/star_particle.h"
#endif
#endif

extern struct star_particle_data
{
  MyIDType PID;

#ifdef METALS
 MyDouble Metals;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
  MyDouble Hsml;
  MyDouble NgbMass;
  MyDouble NgbVolume;
  integertime NgbMinStep;
  int DensityFlag;
  signed char TimeBinStar;
  MyDouble Birthtime;
#ifndef STAR_BY_STAR
  int NumOfStarsInBins[NBINS];
#endif
#endif

#ifdef WINDS 
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
  MyDouble RAD_IonizingHPhotons;
  MyDouble RAD_Ionizing;
#endif

#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
  MyDouble RAD_UVLymanWerner;
  MyDouble RAD_Ultraviolet;
#endif
  
#if defined(RADIATION_PRESSURE)
  MyDouble RAD_Optical;
  MyDouble RAD_Infrared;
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyEject;
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

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
  MyDouble RAD_IonizingRate;
  MyDouble RAD_IonizingLuminosity;
#endif

#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
  MyDouble RAD_UVLymanWernerLuminosity;
  MyDouble RAD_UltravioletLuminosity;
#endif

#if defined(RADIATION_PRESSURE)
  MyDouble RAD_OpticalLuminosity;
  MyDouble RAD_InfraredLuminosity;
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
  int Stage; // 0:preSN, 1:SN, 2:postSN

#ifdef WINDS
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
  MyDouble RAD_IonizingHPhotons;
  MyDouble RAD_Ionizing;
#endif

#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
  MyDouble RAD_UVLymanWerner;
  MyDouble RAD_Ultraviolet;
#endif

#if defined(RADIATION_PRESSURE)
  MyDouble RAD_Optical;
  MyDouble RAD_Infrared;
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