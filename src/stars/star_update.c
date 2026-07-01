#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"


static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}

double gaussian_weight(double r, double h)
{
  double sigma = h / 2.0; 
  double x = r / sigma;
  return exp(-0.5 * x * x);
}

/* Sph loop kernel function */
/*void star_kernel(double u, double hinv3, double hinv4, double *wk, double *dwk)
{
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
}*/

integertime star_timestep(int i)
{ 
 /* Host Hydro Bin */
  double dt_host = (SP[i].HostHydroBin ? (((integertime)1) << SP[i].HostHydroBin) : 0) * All.Timebase_interval;
  
  double dt;
  
  if(dt_host != 0)
    dt = dt_host;
  else 
    dt = TIMEBASE * All.Timebase_interval;

  /* Set a maximum star timestep at 0.01 Myr */
  double dt_star = pow(10,4) / All.cf_UnitTime_in_yr;

  if(dt_star < dt)
    dt = dt_star;

  integertime ti_step = (integertime)(dt / All.Timebase_interval);
  
  return ti_step;
}

void star_update_timesteps(void)
{
  int idx, i;

  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

#if defined(SELFGRAVITY) ||  defined(EXTERNALGRAVITY) || defined(EXACT_GRAVITY_FOR_PARTICLE_TYPE)
      SP[i].TimeBinStar = PPS(i).TimeBinGrav;
#else
      int bin;
      timebins_get_bin_and_do_validity_checks(star_timestep(i), &bin, SP[i].TimeBinStar);
      SP[i].TimeBinStar = bin;
#endif
    }
    
  star_reconstruct_timebins();
  star_update_list_of_active_particles();
}

/* Call this function as the reconstruct_timebins() star version */
void star_reconstruct_timebins(void)
{
  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsStar.TimeBinCount[bin] = 0;
      TimeBinsStar.FirstInTimeBin[bin] = -1;
      TimeBinsStar.LastInTimeBin[bin] = -1;
    }
  
  for(i = 0; i < NumStars; i++)
    {
      
      bin = SP[i].TimeBinStar;
      if(bin >= TIMEBINS)
        continue;

      if(TimeBinsStar.TimeBinCount[bin] > 0)
        {
          TimeBinsStar.PrevInTimeBin[i] = TimeBinsStar.LastInTimeBin[bin];
          TimeBinsStar.NextInTimeBin[i] = -1;
          TimeBinsStar.NextInTimeBin[TimeBinsStar.LastInTimeBin[bin]] = i;
          TimeBinsStar.LastInTimeBin[bin] = i;
        }
      else
        {
          TimeBinsStar.FirstInTimeBin[bin] = TimeBinsStar.LastInTimeBin[bin] = i;
          TimeBinsStar.PrevInTimeBin[i] = TimeBinsStar.NextInTimeBin[i] = -1;
        }
      TimeBinsStar.TimeBinCount[bin]++;
    }
}

/* Call this function after updating the star-timebin */
void star_update_list_of_active_particles(void)
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

  sumup_large_ints(1, &TimeBinsStar.NActiveParticles, &TimeBinsStar.GlobalNActiveParticles);
}

#if defined(WINDS) || defined(SUPERNOVAE)
void feedback_init(struct Mechanical_Feedback_Pack *MFPack)
{
  MFPack->NumEvents = 0;
  MFPack->MaxEvents = 0;
  MFPack->Data = NULL;
}

void feedback_allocate(struct Mechanical_Feedback_Pack *MFPack, int MaxEvents)
{
  MFPack->MaxEvents = MaxEvents;

  MFPack->Data = (Mechanical_Feedback_Data *)
  mymalloc_movable(&MFPack->Data, "Mechanical_Feedback_Events_Pack_Data",
  MFPack->MaxEvents * sizeof(Mechanical_Feedback_Data));
}

void feedback_reallocate(struct Mechanical_Feedback_Pack *MFPack, int NewMaxEvents)
{
  MFPack->MaxEvents = NewMaxEvents;

  MFPack->Data = (Mechanical_Feedback_Data *)
  myrealloc_movable(MFPack->Data,
  MFPack->MaxEvents * sizeof(Mechanical_Feedback_Data));
}
#endif

/* Compute feedback properties of active stars */
void star_prep(void)
{
  TIMER_START(CPU_STARS_PREP);

  int idx, i;
  
  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

      if(SP[i].Active == 0)
        //if(TimeBinSynchronized[SP[i].HostHydroBin])
          {
            SP[i].Active = 1;
            SP[i].PhysicalAge_yr = 0.0;
            SP[i].MassOfStar = PPS(i).Mass;
          }
      
      if(SP[i].Active == 0)
        continue;
      
      MyDouble star_timestep = (SP[i].TimeBinStar ? (((integertime)1) << SP[i].TimeBinStar) : 0) * All.Timebase_interval;

      if(star_timestep == 0)
        terminate("star_timestep == 0!");

      MyDouble star_mass = SP[i].MassOfStar * All.cf_UnitMass_in_Msun;

      /* This sets the timestep of less massive stars at 1 Myr */
      if(star_mass < 8)
        star_timestep = pow(10,6) / All.cf_UnitTime_in_yr;

      /* This deactivates low mass stars */
      if(star_mass < 2)
        star_timestep = TIMEBASE * All.Timebase_interval;

      star_timestep *= All.cf_UnitTime_in_yr;

#ifdef METALS 
      MyDouble star_metallicity = SP[i].Metallicity;
#else 
      MyDouble star_metallicity = 0;
#endif

      SP[i].PhysicalAge_yr += star_timestep;

      Star_Feedback StarFeedback;

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
      StarFeedback = units_for_feedback(star_particle_feedback(i, star_timestep, star_metallicity, SP[i].PhysicalAge_yr));
#elif STAR_PARTICLES == 2     
      StarFeedback = units_for_feedback(star_feedback_compute(star_timestep, star_metallicity, star_mass, SP[i].PhysicalAge_yr));
#endif

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
      SP[i].TimeSN_yr = StarFeedback.TimeSN;
#endif

#if defined(WINDS) || defined(SUPERNOVAE)
      for(int k = 0; k < 3; k++)
        {
          SP[i].WindsAndSN.StarPosition[k] = PPS(i).Pos[k];
          SP[i].WindsAndSN.StarVelocity[k] = PPS(i).Vel[k];
        }
#endif 

#ifdef WINDS
      SP[i].WindsAndSN.MassLoss = StarFeedback.MassLoss;
#ifdef METALS
      SP[i].WindsAndSN.MetalsLoss = StarFeedback.MetalsLoss;
#endif
      SP[i].WindsAndSN.WindMomentum = StarFeedback.WindMomentum;
#endif

#ifdef STAR_RADIATION_ACTIVE      
      for(int w = 0; w < WAVEBANDS; w++)
        {
          SP[i].Radiated[w].Photons = StarFeedback.Radiated[w].Photons;
          SP[i].Radiated[w].Energy = StarFeedback.Radiated[w].Energy;
        }
#endif

#ifdef SUPERNOVAE
      SP[i].WindsAndSN.SN_MassLoss = StarFeedback.SN_MassLoss;
#ifdef METALS
      SP[i].WindsAndSN.SN_MetalsLoss = StarFeedback.SN_MetalsLoss;
#endif
      SP[i].WindsAndSN.SN_EnergyInject = StarFeedback.SN_EnergyInject;
#endif
    }

  TIMER_STOP(CPU_STARS_PREP);
}

void star_perform_end_of_step_physics(void)
{
  int idx, i;

  struct pv_update_data pvd;
  if(All.ComovingIntegrationOn)
    {
      pvd.atime = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

  /* Subtract massloss from stars */
  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

#ifdef WINDS
      PPS(i).Mass -= SP[i].WindsAndSN.MassLoss;
#endif

#ifdef SUPERNOVAE
      PPS(i).Mass -= SP[i].WindsAndSN.SN_MassLoss;
#endif
    }

  /* Dump star injected mass, momentum, and energy into gas */  
  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

#if defined(WINDS) || defined(SUPERNOVAE)
      /* Add mass */ 
      P[i].Mass += SphP[i].StarMassFeed;
      All.StarFeedbackLocal[3] += SphP[i].StarMassFeed;
      
      SphP[i].StarMassFeed = 0;
#ifdef METALS
      /* Add metals */
      SphP[i].GasMetallicity = (SphP[i].GasMetals + SphP[i].StarMetalsFeed) / P[i].Mass;
      sync_conserved_from_primitive(i, METALS_INDEX);
      All.StarFeedbackLocal[4] += SphP[i].StarMetalsFeed;

      SphP[i].StarMetalsFeed = 0;
#endif
#endif
            
#if defined(WINDS) || defined(RADIATION_PRESSURE) || defined(SUPERNOVAE) 
      /* Update momentum */ 
      SphP[i].Momentum[0] += SphP[i].StarMomentumFeed[0];
      SphP[i].Momentum[1] += SphP[i].StarMomentumFeed[1];
      SphP[i].Momentum[2] += SphP[i].StarMomentumFeed[2];
      /* Update velocities */ 
      update_primitive_variables_single(P, SphP, i, &pvd);
      
      /* Set feed flags to zero */
      SphP[i].StarMomentumFeed[0] = SphP[i].StarMomentumFeed[1] = SphP[i].StarMomentumFeed[2] = 0;
#endif

      /* Update total energy */
      SphP[i].Energy += SphP[i].StarEnergyFeed;
      All.StarFeedbackLocal[5] += SphP[i].StarEnergyFeed;
      
      /* Set feed flags to zero */
      SphP[i].StarEnergyFeed = 0;
      
      /* Update internal energy */ 
      update_internal_energy(P, SphP, i, &pvd);
      /* Update pressure */
      set_pressure_of_cell_internal(P, SphP, i);
    } // for(idx...

    MPI_Allreduce(&All.StarFeedbackLocal, &All.StarFeedbackGlobal, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mpi_printf("STARS: Number of active stars = %lld \n", TimeBinsStar.GlobalNActiveParticles);
    mpi_printf("STARS: Mass given by StarParts = %e, Mass taken up by gas particles = %e \n",
    All.StarFeedbackGlobal[0], All.StarFeedbackGlobal[3]);
    mpi_printf("STARS: Metals given by StarParts = %e, Metals taken up by gas particles = %e \n",
    All.StarFeedbackGlobal[1], All.StarFeedbackGlobal[4]);
    mpi_printf("STARS: Energy given by StarParts = %e, Energy taken up by gas particles = %e \n",
    All.StarFeedbackGlobal[2], All.StarFeedbackGlobal[5]);
} 