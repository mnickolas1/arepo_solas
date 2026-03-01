#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../stars/star_tables.h"

#if STAR_PARTICLES == 1
#include "../stars/star_particle.h"
#endif

/*  Feedback tables interpolation */
static inline double linear_interpolation(double x, double x0, double x1, double y0, double y1);
static inline double star_lifetime(int z_idx, double m_val);
static double lifetime(double z_val, double m_val);

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
static double next_SN_time(double z_val, double m_val, double a);
#endif

#if defined(WINDS) || defined(STAR_RADIATION_ACTIVE)
static inline struct star_interpolate interpolate_age(int z_idx, int m_idx, double a);
static struct star_interpolate interpolate_mass(int z_idx, double m_val, double a);
static struct star_interpolate interpolate_metallicity(double z_val, double m_val, double a);
#endif

#ifdef SUPERNOVAE
static inline struct star_interpolate SN_interpolate_mass(int z_idx, double m_val);
static struct star_interpolate SN_interpolate_metallicity(double z_val, double m_val);
#endif

/* Linear interpolation helper function */
static inline double linear_interpolation(double x, double x0, double x1, double y0, double y1) 
{
  // avoid divide by zero
  if(x1 == x0) return y0;

  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

static inline double star_lifetime(int z_idx, double m_val)
{
  // Lifetime for a given metallicity and mass is just the last entry
  // in the corresponding mass-loss table
  
  // Handle mass interpolation
  if(m_val <= M_VALUES[0]) 
    {
      int N_last = N[z_idx][0];
      return Age[z_idx][0][N_last - 1];
    }
    
  if(m_val >= M_VALUES[M_COUNT - 1]) 
    {
      int N_last = N[z_idx][M_COUNT - 1];
      return Age[z_idx][M_COUNT - 1][N_last - 1];
    }

  for(int m = 0; m < M_COUNT - 1; m++) 
    {
      double m0 = M_VALUES[m];
      double m1 = M_VALUES[m + 1];
      if(m_val >= m0 && m_val <= m1) 
        {
          int N0 = N[z_idx][m];
          int N1 = N[z_idx][m + 1];
          double t0 = Age[z_idx][m][N0 - 1];
          double t1 = Age[z_idx][m + 1][N1 - 1];
          return linear_interpolation(m_val, m0, m1, t0, t1);
        }
    }
  terminate("Star_lifetime: failed to bracket mass");
}

static double lifetime(double z_val, double m_val)
{
  if(z_val <= Z_VALUES[0])
    return star_lifetime(0, m_val);
  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return star_lifetime(Z_COUNT - 1, m_val);

  for(int z = 0; z < Z_COUNT - 1; z++) 
    {
      double z0 = Z_VALUES[z];
      double z1 = Z_VALUES[z + 1];
      if(z_val >= z0 && z_val <= z1) 
        {
          double t0 = star_lifetime(z, m_val);
          double t1 = star_lifetime(z + 1, m_val);
          return linear_interpolation(z_val, z0, z1, t0, t1);
        }
    }
  terminate("Lifetime: failed to bracket metallicity");
}

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
static double next_SN_time(double tau, double z_val, double m_val, double a)
{ 
  if(m_val >= 8 && tau > a)
    {
      struct star_interpolate SNfeedback = SN_interpolate_metallicity(z_val, m_val);
      if(SNfeedback.SN_MassLoss > 0.0)
        return tau;
    }
  return MAX_DOUBLE_NUMBER; /* No SN or already past SN */
}
#endif

#if defined(WINDS) || defined(STAR_RADIATION_ACTIVE)
/* Linear interpolation in age */
static inline struct star_interpolate interpolate_age(int z_idx, int m_idx, double a) 
{
  struct star_interpolate feedback = {0};
  
  const double *age = Age[z_idx][m_idx];
  const double *radius = Radius[z_idx][m_idx];
  const double *temperature = Temperature[z_idx][m_idx];
  
#ifdef WINDS
  const double *masslossrate  = MassLossRate[z_idx][m_idx];
#ifdef METALS
  const double *metalslossrate = MetalsLossRate[z_idx][m_idx];
#endif
  const double *windvelocity = WindVelocity[z_idx][m_idx];
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
  const double *RAD_ionizingrate = RAD_IonizingRate[z_idx][m_idx];
  const double *RAD_ionizingluminosity = RAD_IonizingLuminosity[z_idx][m_idx];
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
  const double *RAD_uvlymanwernerluminosity = RAD_UVLymanWernerLuminosity[z_idx][m_idx];
  const double *RAD_ultravioletluminosity = RAD_UltravioletLuminosity[z_idx][m_idx];
#endif
#if defined(RADIATION_PRESSURE)
  const double *RAD_opticalluminosity = RAD_OpticalLuminosity[z_idx][m_idx];
  const double *RAD_infraredluminosity = RAD_InfraredLuminosity[z_idx][m_idx];
#endif

  int n = N[z_idx][m_idx];

  if(a <= age[0])
    {
      feedback.Radius = radius[0];
      feedback.Temperature = temperature[0];

#ifdef WINDS
      feedback.MassLossRate = masslossrate[0];
#ifdef METALS
      feedback.MetalsLossRate = metalslossrate[0];
#endif
      feedback.WindVelocity = windvelocity[0];
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
      feedback.RAD_IonizingRate = RAD_ionizingrate[0];
      feedback.RAD_IonizingLuminosity = RAD_ionizingluminosity[0];
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
      feedback.RAD_UVLymanWernerLuminosity = RAD_uvlymanwernerluminosity[0];
      feedback.RAD_UltravioletLuminosity = RAD_ultravioletluminosity[0];
#endif
#if defined(RADIATION_PRESSURE)
      feedback.RAD_OpticalLuminosity = RAD_opticalluminosity[0];
      feedback.RAD_InfraredLuminosity = RAD_infraredluminosity[0];
#endif
      
      return feedback;
    }       
    
  if(a >= age[n - 1])  
    {
      feedback.Radius = radius[n - 1];
      feedback.Temperature = temperature[n - 1];

#ifdef WINDS
      feedback.MassLossRate = masslossrate[n - 1];
#ifdef METALS
      feedback.MetalsLossRate = metalslossrate[n - 1];
#endif
      feedback.WindVelocity = windvelocity[n - 1];
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
      feedback.RAD_IonizingRate = RAD_ionizingrate[n - 1];
      feedback.RAD_IonizingLuminosity = RAD_ionizingluminosity[n - 1];
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
      feedback.RAD_UVLymanWernerLuminosity = RAD_uvlymanwernerluminosity[n - 1];
      feedback.RAD_UltravioletLuminosity = RAD_ultravioletluminosity[n - 1];
#endif
#if defined(RADIATION_PRESSURE)
      feedback.RAD_OpticalLuminosity = RAD_opticalluminosity[n - 1];
      feedback.RAD_InfraredLuminosity = RAD_infraredluminosity[n - 1];
#endif
      
      return feedback;
    } 
  
  for(int i = 0; i < n - 1; i++)
    {
      if(a >= age[i] && a <= age[i + 1])
        {
          feedback.Radius = linear_interpolation(a, age[i], age[i + 1], radius[i], radius[i + 1]);
          feedback.Temperature = linear_interpolation(a, age[i], age[i + 1], temperature[i], temperature[i + 1]);

#ifdef WINDS
          feedback.MassLossRate = linear_interpolation(a, age[i], age[i + 1], masslossrate[i], masslossrate[i + 1]);
#ifdef METALS
          feedback.MetalsLossRate = linear_interpolation(a, age[i], age[i + 1], metalslossrate[i], metalslossrate[i + 1]);
#endif
          feedback.WindVelocity = linear_interpolation(a, age[i], age[i + 1], windvelocity[i], windvelocity[i + 1]);
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
          feedback.RAD_IonizingRate = linear_interpolation(a, age[i], age[i + 1], RAD_ionizingrate[i], RAD_ionizingrate[i + 1]);
          feedback.RAD_IonizingLuminosity = linear_interpolation(a, age[i], age[i + 1], RAD_ionizingluminosity[i], RAD_ionizingluminosity[i + 1]);
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
          feedback.RAD_UVLymanWernerLuminosity = linear_interpolation(a, age[i], age[i + 1], RAD_uvlymanwernerluminosity[i], RAD_uvlymanwernerluminosity[i + 1]);
          feedback.RAD_UltravioletLuminosity = linear_interpolation(a, age[i], age[i + 1], RAD_ultravioletluminosity[i], RAD_ultravioletluminosity[i + 1]);
#endif
#if defined(RADIATION_PRESSURE)
          feedback.RAD_OpticalLuminosity = linear_interpolation(a, age[i], age[i + 1], RAD_opticalluminosity[i], RAD_opticalluminosity[i + 1]);
          feedback.RAD_InfraredLuminosity = linear_interpolation(a, age[i], age[i + 1], RAD_infraredluminosity[i], RAD_infraredluminosity[i + 1]);
#endif

          return feedback;
        } 
    }
  terminate("Interpolate_age: failed to bracket age");
}

/* Linear interpolation in mass */
static struct star_interpolate interpolate_mass(int z_idx, double m_val, double a) 
{
  if(m_val <= M_VALUES[0])
    return interpolate_age(z_idx, 0, a);
  if(m_val >= M_VALUES[M_COUNT - 1])
    return interpolate_age(z_idx, M_COUNT - 1, a);

  for(int m = 0; m < M_COUNT - 1; m++)
    {
      double m0 = M_VALUES[m];
      double m1 = M_VALUES[m + 1];
      if(m_val >= m0 && m_val <= m1)
        {
          struct star_interpolate feedback0 = interpolate_age(z_idx, m, a);
          struct star_interpolate feedback1 = interpolate_age(z_idx, m + 1, a);
          struct star_interpolate feedback = {0};

          feedback.Radius = linear_interpolation(m_val, m0, m1, feedback0.Radius, feedback1.Radius);
          feedback.Temperature = linear_interpolation(m_val, m0, m1, feedback0.Temperature, feedback1.Temperature);

#ifdef WINDS
          feedback.MassLossRate = linear_interpolation(m_val, m0, m1, feedback0.MassLossRate, feedback1.MassLossRate);
#ifdef METALS
          feedback.MetalsLossRate = linear_interpolation(m_val, m0, m1, feedback0.MetalsLossRate, feedback1.MetalsLossRate);
#endif
          feedback.WindVelocity = linear_interpolation(m_val, m0, m1, feedback0.WindVelocity, feedback1.WindVelocity);
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
          feedback.RAD_IonizingRate = linear_interpolation(m_val, m0, m1, feedback0.RAD_IonizingRate, feedback1.RAD_IonizingRate);
          feedback.RAD_IonizingLuminosity = linear_interpolation(m_val, m0, m1, feedback0.RAD_IonizingLuminosity, feedback1.RAD_IonizingLuminosity);
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
          feedback.RAD_UVLymanWernerLuminosity = linear_interpolation(m_val, m0, m1, feedback0.RAD_UVLymanWernerLuminosity, feedback1.RAD_UVLymanWernerLuminosity);
          feedback.RAD_UltravioletLuminosity = linear_interpolation(m_val, m0, m1, feedback0.RAD_UltravioletLuminosity, feedback1.RAD_UltravioletLuminosity);
#endif
#if defined(RADIATION_PRESSURE)
          feedback.RAD_OpticalLuminosity = linear_interpolation(m_val, m0, m1, feedback0.RAD_OpticalLuminosity, feedback1.RAD_OpticalLuminosity);
          feedback.RAD_InfraredLuminosity = linear_interpolation(m_val, m0, m1, feedback0.RAD_InfraredLuminosity, feedback1.RAD_InfraredLuminosity);
#endif

          return feedback;
        }
    }
  terminate("Interpolate_mass: failed to bracket mass");
}

/* Linear interpolation in metallicity */
static struct star_interpolate interpolate_metallicity(double z_val, double m_val, double a)
{
  if(z_val <= Z_VALUES[0])
    return interpolate_mass(0, m_val, a);
  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return interpolate_mass(Z_COUNT - 1, m_val, a);

  for(int z = 0; z < Z_COUNT - 1; z++)
    {
      double z0 = Z_VALUES[z];
      double z1 = Z_VALUES[z + 1];
      if(z_val >= z0 && z_val <= z1)
        {
          struct star_interpolate feedback0 = interpolate_mass(z, m_val, a);
          struct star_interpolate feedback1 = interpolate_mass(z + 1, m_val, a);
          struct star_interpolate feedback = {0};

          feedback.Radius = linear_interpolation(z_val, z0, z1, feedback0.Radius, feedback1.Radius);
          feedback.Temperature = linear_interpolation(z_val, z0, z1, feedback0.Temperature, feedback1.Temperature);

#ifdef WINDS
          feedback.MassLossRate = linear_interpolation(z_val, z0, z1, feedback0.MassLossRate, feedback1.MassLossRate);
#ifdef METALS
          feedback.MetalsLossRate = linear_interpolation(z_val, z0, z1, feedback0.MetalsLossRate, feedback1.MetalsLossRate);
#endif
          feedback.WindVelocity = linear_interpolation(z_val, z0, z1, feedback0.WindVelocity, feedback1.WindVelocity);
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
          feedback.RAD_IonizingRate = linear_interpolation(z_val, z0, z1, feedback0.RAD_IonizingRate, feedback1.RAD_IonizingRate);
          feedback.RAD_IonizingLuminosity = linear_interpolation(z_val, z0, z1, feedback0.RAD_IonizingLuminosity, feedback1.RAD_IonizingLuminosity);
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
          feedback.RAD_UVLymanWernerLuminosity = linear_interpolation(z_val, z0, z1, feedback0.RAD_UVLymanWernerLuminosity, feedback1.RAD_UVLymanWernerLuminosity);
          feedback.RAD_UltravioletLuminosity = linear_interpolation(z_val, z0, z1, feedback0.RAD_UltravioletLuminosity, feedback1.RAD_UltravioletLuminosity);
#endif
#if defined(RADIATION_PRESSURE)
          feedback.RAD_OpticalLuminosity = linear_interpolation(z_val, z0, z1, feedback0.RAD_OpticalLuminosity, feedback1.RAD_OpticalLuminosity);
          feedback.RAD_InfraredLuminosity = linear_interpolation(z_val, z0, z1, feedback0.RAD_InfraredLuminosity, feedback1.RAD_InfraredLuminosity);
#endif

          return feedback;
        }
    }
  terminate("Interpolate_metallicity: failed to bracket metallicity");
}
#endif

#ifdef SUPERNOVAE
/* Linear interpolation in mass */
static inline struct star_interpolate SN_interpolate_mass(int z_idx, double m_val) 
{
  struct star_interpolate SNfeedback = {0};
  
  const double *SN_massloss  = SN_MassLoss[z_idx];
#ifdef METALS
  const double *SN_metalsloss = SN_MetalsLoss[z_idx];
#endif

  if(m_val <= M_VALUES[0])
    {
      SNfeedback.SN_MassLoss = SN_massloss[0];
#ifdef METALS
      SNfeedback.SN_MetalsLoss = SN_metalsloss[0];
#endif
      SNfeedback.SN_EnergyInject = (SNfeedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
      
      return SNfeedback;
    }       
    
  if(m_val >= M_VALUES[M_COUNT - 1])
    {
      SNfeedback.SN_MassLoss = SN_massloss[M_COUNT - 1];
#ifdef METALS
      SNfeedback.SN_MetalsLoss = SN_metalsloss[M_COUNT - 1];
#endif
      SNfeedback.SN_EnergyInject = (SNfeedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
      
      return SNfeedback;
    } 
  
  for(int m = 0; m < M_COUNT - 1; m++)
    {
      double m0 = M_VALUES[m];
      double m1 = M_VALUES[m + 1];
      if(m_val >= m0 && m_val <= m1)
        {
          SNfeedback.SN_MassLoss = linear_interpolation(m_val, m0, m1, SN_massloss[m], SN_massloss[m + 1]);
#ifdef METALS
          SNfeedback.SN_MetalsLoss = linear_interpolation(m_val, m0, m1, SN_metalsloss[m], SN_metalsloss[m + 1]);
#endif
          SNfeedback.SN_EnergyInject = (SNfeedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
          
          return SNfeedback;
        } 
    }
  terminate("SN_interpolate_mass: failed to bracket mass");
}

/* Linear interpolation in metallicity */
static struct star_interpolate SN_interpolate_metallicity(double z_val, double m_val)
{
  if(z_val <= Z_VALUES[0])
    return SN_interpolate_mass(0, m_val);
  if(z_val >= Z_VALUES[Z_COUNT - 1])
    return SN_interpolate_mass(Z_COUNT - 1, m_val);

  for(int z = 0; z < Z_COUNT - 1; z++)
    {
      double z0 = Z_VALUES[z];
      double z1 = Z_VALUES[z + 1];
      if(z_val >= z0 && z_val <= z1)
        {
          struct star_interpolate SNfeedback0 = SN_interpolate_mass(z, m_val);
          struct star_interpolate SNfeedback1 = SN_interpolate_mass(z + 1, m_val);
          struct star_interpolate SNfeedback = {0};

          SNfeedback.SN_MassLoss = linear_interpolation(z_val, z0, z1, SNfeedback0.SN_MassLoss, SNfeedback1.SN_MassLoss);
#ifdef METALS
          SNfeedback.SN_MetalsLoss = linear_interpolation(z_val, z0, z1, SNfeedback0.SN_MetalsLoss, SNfeedback1.SN_MetalsLoss);
#endif
          SNfeedback.SN_EnergyInject = (SNfeedback.SN_MassLoss > 0.0) ? 1e51 : 0.0;
          
          return SNfeedback;
        }
    }
  terminate("SN_interpolate_metallicity: failed to bracket metallicity");
}
#endif

/* Wrapper function */
struct star_feedback star_feedback_compute(double dt, double z_val, double m_val, double a)
{
  double tau = lifetime(z_val, m_val);
  struct star_feedback star = {0};
  
  star.TimeSN = MAX_DOUBLE_NUMBER;

  if(m_val <= 2)
    return star;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  star.TimeSN = next_SN_time(tau, z_val, m_val, a);
#endif

  if(a < tau)
    {
      star.Stage = 0;

#if defined(WINDS) || defined(STAR_RADIATION_ACTIVE)
      struct star_interpolate feedback = interpolate_metallicity(z_val, m_val, a);

#ifdef WINDS
      star.MassLoss = feedback.MassLossRate * dt;
#ifdef METALS
      star.MetalsLoss = feedback.MetalsLossRate * dt;
#endif
      star.WindMomentum = star.MassLoss * feedback.WindVelocity;
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
      star.RAD_IonizingHPhotons = feedback.RAD_IonizingRate * dt;
      star.RAD_Ionizing = feedback.RAD_IonizingLuminosity * dt;
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
      star.RAD_UVLymanWerner = feedback.RAD_UVLymanWernerLuminosity * dt;
      star.RAD_Ultraviolet = feedback.RAD_UltravioletLuminosity * dt;
#endif
#if defined(RADIATION_PRESSURE)
      star.RAD_Optical = feedback.RAD_OpticalLuminosity * dt;
      star.RAD_Infrared = feedback.RAD_InfraredLuminosity * dt;
#endif
#endif

      return star;
    }

  if(a >= tau && a < (tau+dt)) 
    {
      star.Stage = 1;

#ifdef SUPERNOVAE
      struct star_interpolate SNfeedback = SN_interpolate_metallicity(z_val, m_val);
      star.SN_MassLoss = SNfeedback.SN_MassLoss;
#ifdef METALS
      star.SN_MetalsLoss = SNfeedback.SN_MetalsLoss;
#endif
      star.SN_EnergyInject = SNfeedback.SN_EnergyInject;
#endif

      return star;
    }
  
  if(a >= tau+dt)
    {
      star.Stage = 2;

      return star;
    }
  return star;
}

struct star_feedback units_for_feedback(struct star_feedback star_feedback)
{
#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  star_feedback.TimeSN /= (All.UnitTime_in_s / SEC_PER_YEAR);
#endif

#ifdef WINDS
  star_feedback.MassLoss /= (All.UnitMass_in_g / SOLAR_MASS);
#ifdef METALS
  star_feedback.MetalsLoss /= (All.UnitMass_in_g / SOLAR_MASS);
#endif
  star_feedback.WindMomentum /= ((All.UnitMass_in_g / SOLAR_MASS) * (All.UnitVelocity_in_cm_per_s / 1.e5));
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
  star_feedback.RAD_Ionizing /= (All.UnitEnergy_in_cgs);
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
  star_feedback.RAD_UVLymanWerner /= (All.UnitEnergy_in_cgs);
  star_feedback.RAD_Ultraviolet /= (All.UnitEnergy_in_cgs);
#endif
#if defined(RADIATION_PRESSURE)
  star_feedback.RAD_Optical /= (All.UnitEnergy_in_cgs);
  star_feedback.RAD_Infrared /= (All.UnitEnergy_in_cgs);
#endif

#ifdef SUPERNOVAE
  star_feedback.SN_MassLoss /= (All.UnitMass_in_g / SOLAR_MASS);
#ifdef METALS
  star_feedback.SN_MetalsLoss /= (All.UnitMass_in_g / SOLAR_MASS);
#endif
  star_feedback.SN_EnergyInject /= (All.UnitEnergy_in_cgs);
#endif

  return star_feedback;
}

#if STAR_PARTICLES == 1
struct star_feedback star_particle_feedback(int index, double dt, double z, double a)
{  
  int i, Nstars;
  double m;
  struct star_feedback star_particle = {0};
  
  star_particle.TimeSN = MAX_DOUBLE_NUMBER; 

  // Add feedback contributions for each bin 
  for(i = 0; i < NBINS; i++) 
    {
      if(SP[index].NumOfStarsInBins[i] == 0)
        continue;

      Nstars = SP[index].NumOfStarsInBins[i];
      
      m = StarMeanMassInBins[i]; 

      struct star_feedback star = star_feedback_compute(dt, z, m, a);

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
      if(star.TimeSN < star_particle.TimeSN)
      star_particle.TimeSN = star.TimeSN;
#endif

      switch(star.Stage)
        {
          case 0:

#ifdef WINDS
          star_particle.MassLoss += Nstars * star.MassLoss;
#ifdef METALS
          star_particle.MetalsLoss += Nstars * star.MetalsLoss;
#endif
          star_particle.WindMomentum += Nstars * star.WindMomentum;
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
          star_particle.RAD_IonizingHPhotons += Nstars * star.RAD_IonizingHPhotons;
          star_particle.RAD_Ionizing += Nstars * star.RAD_Ionizing;
#endif
#if defined(PHOTOELECTRIC_HEATING) || defined(RADIATION_PRESSURE)
          star_particle.RAD_UVLymanWerner += Nstars * star.RAD_UVLymanWerner;
          star_particle.RAD_Ultraviolet += Nstars * star.RAD_Ultraviolet;
#endif
#if defined(RADIATION_PRESSURE)
          star_particle.RAD_Optical += Nstars * star.RAD_Optical;
          star_particle.RAD_Infrared += Nstars * star.RAD_Infrared;
#endif

          break;

          case 1:

#ifdef SUPERNOVAE
          star_particle.SN_MassLoss += Nstars * star.SN_MassLoss;
#ifdef METALS
          star_particle.SN_MetalsLoss += Nstars * star.SN_MetalsLoss;
#endif
          star_particle.SN_EnergyInject += Nstars * star.SN_EnergyInject;   
#endif

          break;

          case 2:

          break;
        }
    }
  return star_particle;
}
#endif