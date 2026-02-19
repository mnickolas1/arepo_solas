#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"

static int int_compare(const void *a, const void *b);

/* Sph loop kernel function */
static void star_kernel(double u, double hinv3, double hinv4, double *wk, double *dwk)
{
  // Cubic spline
  double K_norm = 8.0 / M_PI;

  if(u < 0.5)
    {
      *dwk = u * (18.0 * u - 12.0);
      
      *wk = (1.0 + 6.0 * (u - 1.0) * u * u);
    }
  else
    {
      double t1 = (1.0 - u);
      double t2 = t1 * t1;
      
      *dwk = -6.0 * t2;
      
      *wk = 2.0 * t2 * t1;
    }

  *dwk *= K_norm * hinv4;
  
  *wk  *= K_norm * hinv3;
}

/* Compute feedback properties of active stars */
void star_prep(void)
{
  int i, idx;
  
  for(idx=0; idx<TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];
      
      MyDouble star_timestep = All.TimeStep;
      MyDouble star_age = All.Time - SP[i].Birthtime;
      MyDouble star_mass = PPS(i).Mass; 
      MyDouble star_metals = SP[i].Metals;
      MyDouble star_metallicity = star_metals/star_mass;
    
      // Convert units (-> solar masses and years) (cosmological?)
      star_timestep *= (All.cf_atime / All.cf_time_hubble_a) * (All.UnitTime_in_s / SEC_PER_YEAR);
      star_age *= (All.UnitTime_in_s / SEC_PER_YEAR);
      star_mass *= (All.UnitMass_in_g / SOLAR_MASS);


      struct star_feedback star_feedback;

#ifndef STAR_BY_STAR
      star_feedback = units_for_feedback(star_particle_feedback(star_timestep, star_metallicity, star_mass, star_age));
#else     
      star_feedback = units_for_feedback(star_feedback_compute(star_timestep, star_metallicity, star_mass, star_age));
#endif

#ifdef WINDS
      SP[i].MassLoss = star_feedback.MassLoss;
#ifdef METALS
      SP[i].MetalsLoss = star_feedback.MetalsLoss;
#endif
      SP[i].WindMomentum = star_feedback.WindMomentum;
#endif

#ifdef SUPERNOVAE
      SP[i].SN_MassLoss = star_feedback.SN_MassLoss;
#ifdef METALS
      SP[i].SN_MetalsLoss = star_feedback.SN_MetalsLoss;
#endif
      SP[i].SN_EnergyEject = star_feedback.SN_EnergyInject;
#endif

#if defined(PHOTOIONIZATION) || defined(RADIATION_PRESSURE)
      SP[i].RAD_IonizingHPhotons = star_feedback.RAD_IonizingHPhotons;
      SP[i].RAD_Ionizing = star_feedback.RAD_Ionizing;
#endif
#if defined(PHOTOELECTRIC) || defined(RADIATION_PRESSURE)
      SP[i].RAD_UVLymanWerner = star_feedback.RAD_UVLymanWerner;
      SP[i].RAD_Ultraviolet = star_feedback.RAD_Ultraviolet;
#endif
#if defined(RADIATION_PRESSURE)
      SP[i].RAD_Optical = star_feedback.RAD_Optical;
      SP[i].RAD_Infrared = star_feedback.RAD_Infrared;
#endif
    }
}

/* Get timestep for star based on smallest between ngbs */
integertime get_timestep_star(int p)
{ 
  // Star particles always active for winds
  return 0;
  
  // return (integertime)(1e-2 / All.Timebase_interval);
}

void update_star_timesteps(void)
{
  int i, bin;
  integertime ti_step;

  for(i = 0; i < NumStars; i++)
    { 
      ti_step = get_timestep_star(i);
    
      bin = get_timestep_bin(ti_step);

      SP[i].TimeBinStar = bin;
    }
  reconstruct_star_timebins();
  update_list_of_active_star_particles();
}

/* Call this function as the reconstruct_timebins() star version */
void reconstruct_star_timebins(void)
{
  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsStar.TimeBinCount[bin]   = 0;
      TimeBinsStar.FirstInTimeBin[bin] = -1;
      TimeBinsStar.LastInTimeBin[bin]  = -1;
    }
  
  for(i = 0; i < NumStars; i++)
    {
      
      bin = SP[i].TimeBinStar;
      if(bin >= TIMEBINS)
        continue;

      if(TimeBinsStar.TimeBinCount[bin] > 0)
        {
          TimeBinsStar.PrevInTimeBin[i]                                  = TimeBinsStar.LastInTimeBin[bin];
          TimeBinsStar.NextInTimeBin[i]                                  = -1;
          TimeBinsStar.NextInTimeBin[TimeBinsStar.LastInTimeBin[bin]]    = i;
          TimeBinsStar.LastInTimeBin[bin]                                = i;
        }
      else
        {
          TimeBinsStar.FirstInTimeBin[bin] = TimeBinsStar.LastInTimeBin[bin] = i;
          TimeBinsStar.PrevInTimeBin[i] = TimeBinsStar.NextInTimeBin[i] = -1;
        }
      TimeBinsStar.TimeBinCount[bin]++;
    }
}

/* Call this function after updating the star-timebin to the ngb condition */
void update_list_of_active_star_particles(void)
{
  int i, n;
  TimeBinsStar.NActiveParticles = 0;
  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n]) 
        {
          for(i = TimeBinsStar.FirstInTimeBin[n]; i >= 0; i = TimeBinsStar.NextInTimeBin[i])
            {
              TimeBinsStar.ActiveParticleList[TimeBinsStar.NActiveParticles] = i;
              TimeBinsStar.NActiveParticles++;  
            }
        }
    }

    mysort(TimeBinsStar.ActiveParticleList, TimeBinsStar.NActiveParticles, sizeof(int), int_compare);
}

void perform_end_of_step_star_physics(void)
{
  int idx, i;
  double pj, p0;
  double kick_vector[3];
    
  struct pv_update_data pvd;
  if(All.ComovingIntegrationOn)
    {
      pvd.atime    = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;
    
  // Inject feedback to ngb cells
  if(All.Time >= All.FeedbackTime)
    {        
      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
            continue;
            
          // Dump mass, momentum and energy injected by stars 
#if defined(WINDS) || defined(SUPERNOVAE)
          // Add mass 
          P[i].Mass += SphP[i].MassFeed;
#ifdef METALS
          // Add metals
          SphP[i].Metals += SphP[i].MetalsFeed;
#endif
#endif
            
#if defined(WINDS) || defined(SUPERNOVAE) || defined (RADIATION_PRESSURE)
          // Calculate kick
          kick_vector[0] = SphP[i].MomentumKickVector[0];
          kick_vector[1] = SphP[i].MomentumKickVector[1];
          kick_vector[2] = SphP[i].MomentumKickVector[2];

#ifdef WINDS
          if(SphP[i].MomentumFeed > 0)
            {
              // Momentum conserving wind 
              pj = SphP[i].MomentumFeed;
                
              // Update momentum 
              SphP[i].Momentum[0] += kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              SphP[i].Momentum[1] += kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              SphP[i].Momentum[2] += kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
                
              All.EnergyExchange[3] += SphP[i].MomentumFeed;
                
              // Update velocities 
              update_primitive_variables_single(P, SphP, i, &pvd);
              // Update total energy 
              SphP[i].Energy = SphP[i].Utherm * P[i].Mass +
                0.5 * P[i].Mass * (pow(P[i].Vel[0], 2) + pow(P[i].Vel[1], 2) + pow(P[i].Vel[2], 2));
              // Update internal energy 
              update_internal_energy(P, SphP, i, &pvd);
              // Update pressure
              set_pressure_of_cell_internal(P, SphP, i);
              // Set feed flags to zero
              SphP[i].MomentumFeed = 0;
            }
#endif

#ifdef SUPERNOVAE
          if(SphP[i].ThermalEnergyFeed || SphP[i].KineticEnergyFeed> 0)
            {
              // Energy conserving supernova
             
              // Add kinetic energy 
              SphP[i].Energy += SphP[i].KineticEnergyFeed * All.cf_atime*All.cf_atime;
                
              // Calculate momentum feed exactly so energy is conserved 
              //-> we need to do this here so that particle properties don't change between loading the buffer and emptying it
              p0 = sqrt(pow(SphP[i].Momentum[0], 2) + pow(SphP[i].Momentum[1], 2) + pow(SphP[i].Momentum[2], 2));
                
              pj = sqrt(2 * P[i].Mass * (SphP[i].Energy - (P[i].Mass/*-SphP[i].MassFeed*/)*SphP[i].Utherm*All.cf_atime*All.cf_atime)) - p0;
  
              // Update total energy 
              SphP[i].Energy += SphP[i].ThermalEnergyFeed * All.cf_atime*All.cf_atime;
              All.EnergyExchange[5] += SphP[i].ThermalEnergyFeed + SphP[i].KineticEnergyFeed;
              // Update momentum 
              SphP[i].Momentum[0] += kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              SphP[i].Momentum[1] += kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              SphP[i].Momentum[2] += kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              // Update velocities 
              update_primitive_variables_single(P, SphP, i, &pvd);
              // Update internal energy 
              update_internal_energy(P, SphP, i, &pvd);
              // Update pressure 
              set_pressure_of_cell_internal(P, SphP, i);
              // Set feed flags to zero
              SphP[i].ThermalEnergyFeed = SphP[i].KineticEnergyFeed = 0;
            }
#endif
#endif
        } // for(idx...
        
    } // if(All.Time >= All.FeedbackTime)

    MPI_Allreduce(&All.EnergyExchange, &All.EnergyExchangeTot, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  
    mpi_printf("STARS: Momentum given by StarParts = %e, Momentum taken up by gas particles = %e \n",
               All.EnergyExchangeTot[2] * All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s,
               All.EnergyExchangeTot[3] * All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s);
    mpi_printf("STARS: Energy given by StarParts = %e, Energy taken up by gas particles = %e \n",
               All.EnergyExchangeTot[4] * All.UnitEnergy_in_cgs, All.EnergyExchangeTot[5] * All.UnitEnergy_in_cgs);
    
} 

static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}