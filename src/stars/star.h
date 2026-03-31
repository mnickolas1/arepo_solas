#ifndef STAR_H
#define STAR_H

#include "../main/allvars.h"

#define ALLOC_STAR_ROOM 64
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
|| defined(PHOTOIONIZATION) || defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE) \
|| defined(SUPERNOVAE)
#define STAR_FEEDBACK_ACTIVE
#endif

#if defined(PHOTOIONIZATION) || defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
#define STAR_RADIATION_ACTIVE
#endif

#if defined(STAR_FEEDBACK_ACTIVE) 
#ifndef STAR_PARTICLES
#error "We cannot run star feedback simulations without a star particles model!"
#endif
#endif

#if defined(STAR_PARTICLES) || defined(STAR_FEEDBACK_ACTIVE)
#ifndef STARS
#error "We cannot run star feedback simulations without stars!"
#endif
#endif

#if defined(STAR_PARTICLES) 
#if !defined(USE_SFR) && !defined(STAR_FEEDBACK_ACTIVE) 
#error "We need star formation/feedback for star particles!"
#endif
#endif

#if defined(STAR_RADIATION_ACTIVE)
#ifndef USE_GRACKLE
#error "We cannot run star radiation simulations without GRACKLE!"
#endif
#endif

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
  MyDouble Metals;
#endif

#if STAR_PARTICLES < 2
  int NumOfStarsInBins[NBINS];
#endif

#if STAR_PARTICLES == 2
  MyDouble MassOfStar;
#endif

#ifdef INDIVIDUAL_STAR_BY_STAR_FORMATION
  MyDouble MassToDrain;
#endif

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeSN;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
  int Active;
  MyDouble Hsml;
  MyDouble NgbMass;
  MyDouble NgbVolume;
  int NgbMaxBin;
  int DensityFlag;
  signed char TimeBinStar;
  MyDouble Birthtime;
#endif

#ifdef WINDS 
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE
  MyDouble LUM[WAVEBANDS];
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