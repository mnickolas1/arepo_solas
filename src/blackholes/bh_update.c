#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "../main/allvars.h"
#include "../main/proto.h"

static int int_compare(const void *a, const void *b);

/* Sph loop kernel function */
void bh_kernel(double u, double hinv3, double hinv4, double *wk, double *dwk)
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

#ifdef BLACKHOLE_ACCRETION_ACTIVE
#ifdef BONDI_ACCRETION 
/* Calculate bondi accretion rate */
void update_bh_accretion_rate(void)
{
  
  int i;
  double density, pressure, sound_speed, velocity_gas_norm;
  double denominator, denominator_inv, BondiRate, EddingtonRate;
  double accretion_rate, acc_rate_for_print;

  accretion_rate = acc_rate_for_print = 0;

  for(i = 0; i < NumBhs; i++)
    {
      // Get pressure
      if(BhP[i].Density>0)
        {  
          density = BhP[i].Density;
          pressure = GAMMA_MINUS1 * density * BhP[i].InternalEnergyGas;

          // Get soundspeed
          sound_speed = sqrt(GAMMA * pressure / density);
      
          velocity_gas_norm = sqrt(BhP[i].VelocityGas[0]*BhP[i].VelocityGas[0] + 
            BhP[i].VelocityGas[1]*BhP[i].VelocityGas[1] + BhP[i].VelocityGas[2]*BhP[i].VelocityGas[2]); 

          denominator = (sound_speed*sound_speed + velocity_gas_norm*velocity_gas_norm);
          if(denominator > 0)
            {
              denominator_inv = 1. / sqrt(denominator);
              BondiRate = 4. * M_PI * All.G * All.G * PPB(i).Mass * PPB(i).Mass * density *
              denominator_inv * denominator_inv * denominator_inv;
            }
          else
            terminate("Invalid denominator in Bondi Accretion Rate");
        }
      else
        BondiRate = 0;
  
      // Limit by Eddington accretion rate 
      EddingtonRate = 4. * M_PI * GRAVITY * (PPB(i).Mass * All.UnitMass_in_g) * PROTONMASS / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *=  (All.UnitTime_in_s / All.UnitMass_in_g);
      accretion_rate = fmin(BondiRate, EddingtonRate);
      
      BhP[i].AccretionRate  = accretion_rate;
    }
 
  MPI_Allreduce(&accretion_rate, &acc_rate_for_print, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  mpi_printf("BLACK_HOLES: Black hole accretion rate: %e \n", acc_rate_for_print);
}
#endif
#endif

/* Get timestep for bh based on smallest between ngbs */
integertime get_timestep_bh(int p)
{ 
  return BhP[p].NgbMinStep;
}

/* Update bh-timestep at prior_mesh_construction */
void update_bh_timesteps(void)
{
  int i, bin;
  integertime ti_step;

  for(i = 0; i < NumBhs; i++)
    { 
      ti_step = get_timestep_bh(i);
      //binold = BhP[i].TimeBinBh;
      bin = get_timestep_bin(ti_step);

      //timebin_move_particle(&TimeBinsBh, i, binold, bin);
      BhP[i].TimeBinBh = bin;
    }
  reconstruct_bh_timebins();
  update_list_of_active_bh_particles();
}

/* Call this function as the reconstruct_timebins() bh version */
void reconstruct_bh_timebins(void)
{
  int i, bin;

  for(bin = 0; bin < TIMEBINS; bin++)
    {
      TimeBinsBh.TimeBinCount[bin]   = 0;
      TimeBinsBh.FirstInTimeBin[bin] = -1;
      TimeBinsBh.LastInTimeBin[bin]  = -1;
    }
  
  for(i = 0; i < NumBhs; i++)
    {
      
      bin = BhP[i].TimeBinBh;
      if(bin >= TIMEBINS)
        continue;

      if(TimeBinsBh.TimeBinCount[bin] > 0)
        {
          TimeBinsBh.PrevInTimeBin[i]                                  = TimeBinsBh.LastInTimeBin[bin];
          TimeBinsBh.NextInTimeBin[i]                                  = -1;
          TimeBinsBh.NextInTimeBin[TimeBinsBh.LastInTimeBin[bin]]      = i;
          TimeBinsBh.LastInTimeBin[bin]                                = i;
        }
      else
        {
          TimeBinsBh.FirstInTimeBin[bin] = TimeBinsBh.LastInTimeBin[bin] = i;
          TimeBinsBh.PrevInTimeBin[i] = TimeBinsBh.NextInTimeBin[i] = -1;
        }
      TimeBinsBh.TimeBinCount[bin]++;
    }
}

/* Call this function after updating the bh-timebin to the ngb condition */
void update_list_of_active_bh_particles(void)
{
  int i, n;
  TimeBinsBh.NActiveParticles = 0;
  for(n = 0; n < TIMEBINS; n++)
    {
      if(TimeBinSynchronized[n]) 
        {
          for(i = TimeBinsBh.FirstInTimeBin[n]; i >= 0; i = TimeBinsBh.NextInTimeBin[i])
            {
              TimeBinsBh.ActiveParticleList[TimeBinsBh.NActiveParticles] = i;
              TimeBinsBh.NActiveParticles++;  
            }
        }
    }

    mysort(TimeBinsBh.ActiveParticleList, TimeBinsBh.NActiveParticles, sizeof(int), int_compare);

  /*n = 1;
  int in;
  long long out;

  in = TimeBinsBh.NActiveParticles;

  sumup_large_ints(n, &in, &out);

  TimeBinsBh.GlobalNActiveParticles = out;*/
}

void perform_end_of_step_bh_physics(void)
{
  int idx, i;
  double pj, p0;
  double kick_vector[3];
    
#ifdef BLACKHOLE_ACCRETION_ACTIVE
#ifdef BONDI_ACCRETION
  int j, bin;
  double dt;
  // Accrete mass, angular momentum onto the bh and drain ngb cells
  for(i=0; i<NumBhs; i++)
    {
      bin = BhP[i].TimeBinBh;
      dt  = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;
      PPB(i).Mass += (1-All.Epsilon_r) * BhP[i].AccretionRate * dt;
      BhP[i].AngularMomentum[0] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[0];
      BhP[i].AngularMomentum[1] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[1];
      BhP[i].AngularMomentum[2] += BhP[i].AccretionRate * dt * BhP[i].VelocityGasCircular[2];
      for(j=0; j<NumGas; j++)
        {
          if(SphP[j].MassDrain > 0)
            {
              if(P[j].Mass - SphP[j].MassDrain < 0.1*P[j].Mass)
                {
                  P[j].Mass -= 0.9*P[j].Mass;
                  BhP[i].MassToDrain += SphP[j].MassDrain - 0.9*P[j].Mass;
                  // We're also losing thermal and kinetic energy & momentum 
                    
                  // Update total energy 
                  SphP[j].Energy *= 0.1;
                    
                  // Update momentum 
                  SphP[j].Momentum[0] *= 0.1;
                  SphP[j].Momentum[1] *= 0.1;
                  SphP[j].Momentum[2] *= 0.1;
                }
              else
                {
                  P[j].Mass -= SphP[j].MassDrain;
                    
                  // Update total energy 
                  SphP[j].Energy *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);
                    
                  // Update momentum 
                  SphP[j].Momentum[0] *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);
                  SphP[j].Momentum[1] *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);
                  SphP[j].Momentum[2] *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);
                }
              SphP[j].MassDrain = 0;
            }
        }
    }
#endif
    
#ifdef INFALL_ACCRETION
  for(i=0; i<NumBhs; i++)
    {
      PPB(i).Mass += (1-All.Epsilon_r) * BhP[i].Accretion;
      BhP[i].Accretion = 0;
    }
#endif
#endif

#ifdef BLACKHOLE_FEEDBACK_ACTIVE   
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
      if(All.FeedbackFlag > 0)
        {
          for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              i = TimeBinsHydro.ActiveParticleList[idx];
              if(i < 0)
                continue;
                
              // Dump mass, momentum and energy injected by bh 
              if(SphP[i].ThermalFeed > 0 || SphP[i].KineticFeed > 0)
                {
                  // Add mass
                  P[i].Mass += SphP[i].MassLoading;
                    
                  // Add kinetic energy */
                  SphP[i].Energy += SphP[i].KineticFeed * All.cf_atime*All.cf_atime;
                    
                  // Calculate momentum feed exactly so energy is conserved */
                  //-> we need to do this here so that particle properties don't change between loading the buffer and emptying it*/
                  kick_vector[0] = SphP[i].BhKickVector[0];
                  kick_vector[1] = SphP[i].BhKickVector[1];
                  kick_vector[2] = SphP[i].BhKickVector[2];
                    
                  p0 = sqrt(pow(SphP[i].Momentum[0], 2) + pow(SphP[i].Momentum[1], 2) + pow(SphP[i].Momentum[2], 2));
                    
                  pj = sqrt(2 * P[i].Mass * (SphP[i].Energy - (P[i].Mass-SphP[i].MassLoading)*SphP[i].Utherm*All.cf_atime*All.cf_atime)) - p0;

                  // Update total energy 
                  SphP[i].Energy += SphP[i].ThermalFeed * All.cf_atime*All.cf_atime;
                  All.EnergyExchange[1] += SphP[i].ThermalFeed + SphP[i].KineticFeed;
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
                  SphP[i].ThermalFeed = SphP[i].KineticFeed = SphP[i].MassLoading = 0;
                  momentum_kick[0] = momentum_kick[1] = momentum_kick[2] = 0;
#ifdef PASSIVE_SCALARS
                  // Tracer field advected passively 
                  SphP[i].PScalars[0] = 1;
                  SphP[i].PConservedScalars[0] = P[i].Mass;
#endif
                }
            }
#ifdef BURST_MODE
            All.FeedbackFlag = -1;
#endif
        } // if(All.FeedbackFlag>0)
    } // if(All.Time >= All.FeedbackTime)
        
    MPI_Allreduce(&All.EnergyExchange, &All.EnergyExchangeTot, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); 
  
    mpi_printf("BLACK_HOLES: Energy given by BH = %e, Energy taken up by gas particles = %e \n",
               All.EnergyExchangeTot[0] * All.UnitEnergy_in_cgs, All.EnergyExchangeTot[1] * All.UnitEnergy_in_cgs);
#endif /* ifdef BLACKHOLES */
    
#ifdef BURST_MODE
    if(All.EnergyExchangeTot[0] - All.EnergyExchangeTot[1] > 10)
        All.FeedbackFlag = 1;
#endif
} // perform_end_of_step_physics(void)

static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}