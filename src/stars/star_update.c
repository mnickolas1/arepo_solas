#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"


static int int_compare(const void *a, const void *b);

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


int get_timestep_star(int p)
{ 
  double dt_grav = (PPS(p).TimeBinGrav ? (((integertime)1) << PPS(p).TimeBinGrav) : 0) * All.Timebase_interval;
  double dt_ngbmax = (SP[p].NgbMaxBin ? (((integertime)1) << SP[p].NgbMaxBin) : 0) * All.Timebase_interval;
  double dt_star = pow(10,4) * SEC_PER_YEAR / All.UnitTime_in_s;

  double dt = dt_grav;

  if(dt_ngbmax < dt)
    dt = dt_ngbmax;

  if(dt_star < dt)
    dt = dt_star;

  dt *= All.cf_hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

  if(dt < All.MinSizeTimestep)
    {
#ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP
      dt = All.MinSizeTimestep;
#else  /* #ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP */
      print_particle_info(p);
      terminate("Timestep dt=%g below All.MinSizeTimestep=%g", dt, All.MinSizeTimestep);
#endif /* #ifdef NOSTOP_WHEN_BELOW_MINTIMESTEP #else */
    }

  integertime ti_step = (integertime)(dt / All.Timebase_interval);

  validate_timestep(dt, ti_step, p);
  
  return ti_step;
}

void update_star_timesteps(void)
{
  int idx, i;
  integertime ti_step;

for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

      ti_step = get_timestep_star(i);
    
      SP[i].TimeBinStar = get_timestep_bin(ti_step);
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

/* Compute feedback properties of active stars */
void star_prep(void)
{
  int idx, i;
  
  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];

      if(SP[i].Active == 0)
        if(TimeBinSynchronized[SP[i].NgbMaxBin])
          {
            SP[i].Active = 1;
            
            // need to reset the age to when star first starts injecting feedback
            SP[i].Birthtime = All.Time; 
          }
      
      if(SP[i].Active == 0)
        continue;
      
      MyDouble star_timestep = (SP[i].TimeBinStar ? (((integertime)1) << SP[i].TimeBinStar) : 0) * All.Timebase_interval;
      MyDouble star_age = All.Time - SP[i].Birthtime;
      MyDouble star_mass = PPS(i).Mass;
#ifdef METALS 
      MyDouble star_metals = SP[i].Metals;
#else 
      MyDouble star_metals = 0;
#endif
      MyDouble star_metallicity = star_metals/star_mass;
    
      // Convert units (-> solar masses and years)
      star_timestep *= (All.cf_atime / All.cf_time_hubble_a) * (All.UnitTime_in_s / SEC_PER_YEAR);
      star_age *= (All.cf_atime / All.cf_time_hubble_a) * (All.UnitTime_in_s / SEC_PER_YEAR);
      star_mass *= (All.UnitMass_in_g / SOLAR_MASS);

      struct star_feedback star_feedback;

#if STAR_PARTICLES == 1
      star_feedback = units_for_feedback(star_particle_feedback(i, star_timestep, star_metallicity, star_age));
#elif STAR_PARTICLES == 2     
      star_feedback = units_for_feedback(star_feedback_compute(star_timestep, star_metallicity, star_mass, star_age));
#endif

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
     SP[i].TimeSN = star_feedback.TimeSN;
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
      SP[i].SN_EnergyInject = star_feedback.SN_EnergyInject;
#endif

#ifdef STAR_RADIATION_ACTIVE
      double LUM[WAVEBANDS];
      LUM[0] = star_feedback.RAD_IonizingHPhotons;
      LUM[1] = star_feedback.RAD_Ionizing;
      LUM[2] = star_feedback.RAD_UVLymanWerner;
      LUM[3] = star_feedback.RAD_Ultraviolet;
      LUM[4] = star_feedback.RAD_Optical;
      LUM[5] = star_feedback.RAD_Infrared;
      
      for(int w = 0; w < WAVEBANDS; w++)
        SP[i].LUM[w] = LUM[w];
#endif
    }
}

void perform_end_of_step_star_physics(void)
{
  int idx, i;

  struct pv_update_data pvd;
  if(All.ComovingIntegrationOn)
    {
      pvd.atime    = All.Time;
      pvd.hubble_a = hubble_function(All.Time);
      pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
    }
  else
    pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

/*#if defined(WINDS) || defined(SUPERNOVAE)
  for(idx=0, npleft=0; idx<TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];
    }
#endif*/

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
          double Mold = P[i].Mass;

          P[i].Mass += SphP[i].StarMassFeed;
          All.StarFeedbackLocal[4] += SphP[i].StarMassFeed;
          SphP[i].StarMassFeed = 0;
#ifdef METALS
          // Add metals
          SphP[i].GasMetallicity = (Mold * SphP[i].GasMetallicity + SphP[i].StarMetalsFeed)/(P[i].Mass);
          All.StarFeedbackLocal[5] += SphP[i].StarMetalsFeed;
          SphP[i].StarMetalsFeed = 0;
#endif
#endif
            
#if defined(WINDS) || defined(RADIATION_PRESSURE) || defined(SUPERNOVAE) 
          // Update momentum 
          SphP[i].Momentum[0] += SphP[i].StarMomentumFeed[0] * All.cf_atime;
          SphP[i].Momentum[1] += SphP[i].StarMomentumFeed[1] * All.cf_atime;
          SphP[i].Momentum[2] += SphP[i].StarMomentumFeed[2] * All.cf_atime;
          All.StarFeedbackLocal[6] += sqrt(pow(SphP[i].StarMomentumFeed[0],2) + 
          pow(SphP[i].StarMomentumFeed[1],2) + pow(SphP[i].StarMomentumFeed[2],2));   
          // Update velocities 
          update_primitive_variables_single(P, SphP, i, &pvd);
          // Set feed flags to zero
          SphP[i].StarMomentumFeed[0] = SphP[i].StarMomentumFeed[1] = SphP[i].StarMomentumFeed[2] = 0;
#endif

#if defined(RADIATION_PRESSURE) && !defined(WINDS) && !defined(PHOTOIONIZATION) && !defined(PHOTOELECTRIC_HEATING) && !defined(SUPERNOVAE) 
          // Update total energy
          double Eold = SphP[i].Energy;
          SphP[i].Energy = 0.5*P[i].Mass*(P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2])
          + P[i].Mass * SphP[i].Utherm;
          All.StarFeedbackLocal[7] += SphP[i].Energy - Eold;
#else if       
          // Update total energy
          SphP[i].Energy += SphP[i].StarEnergyFeed * All.cf_atime * All.cf_atime;
          All.StarFeedbackLocal[7] += SphP[i].StarEnergyFeed;
          // Set feed flags to zero
          SphP[i].StarEnergyFeed = 0;
#endif
          // Update internal energy 
          update_internal_energy(P, SphP, i, &pvd);
          // Update pressure
          set_pressure_of_cell_internal(P, SphP, i);

        } // for(idx...
        
    } // if(All.Time >= All.FeedbackTime)

    MPI_Allreduce(&All.StarFeedbackLocal, &All.StarFeedbackGlobal, 8, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mpi_printf("STARS: Number of StarParts = %d \n", All.TotNumStars);
    mpi_printf("STARS: Mass given by StarParts = %e, Mass taken up by gas particles = %e \n",
               All.StarFeedbackGlobal[0], All.StarFeedbackGlobal[4]);
    mpi_printf("STARS: Metals given by StarParts = %e, Metals taken up by gas particles = %e \n",
               All.StarFeedbackGlobal[1], All.StarFeedbackGlobal[5]);
    //mpi_printf("STARS: Momentum given by StarParts = %e, Momentum taken up by gas particles = %e \n",
               //All.StarFeedbackGlobal[2] * All.cf_atime, All.StarFeedbackGlobal[6] * All.cf_atime);
    mpi_printf("STARS: Energy given by StarParts = %e, Energy taken up by gas particles = %e \n",
               All.StarFeedbackGlobal[3] * All.cf_atime * All.cf_atime, All.StarFeedbackGlobal[7] * All.cf_atime * All.cf_atime);   
} 

static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}