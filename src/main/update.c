#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../../celib/src/config.h"

#if defined(STARS) || defined(BLACKHOLES)
/* THIS PART ADAPTED FROM GADGET4 */
/* fall back to cubic spline kernel */
#if !defined(CUBIC_SPLINE_KERNEL) && !defined(WENDLAND_C2_KERNEL) && !defined(WENDLAND_C4_KERNEL) && !defined(WENDLAND_C6_KERNEL)
#define CUBIC_SPLINE_KERNEL 
#endif

/* fall back to three dimensions */
#if !defined(TWODIMS) && !defined(ONEDIMS)
#define THREEDIMS
#endif

/* Norms */
#ifdef CUBIC_SPLINE_KERNEL

#ifdef THREEDIMS
#define K_NORM (8.0 / M_PI) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define K_NORM (40.0 / (7.0 * M_PI)) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define K_NORM (4.0 / 3.0) /*!< For 1D-normalized kernel */
#endif

#endif /* CUBIC_SPLINE_KERNEL */

#ifdef WENDLAND_C2_KERNEL

#ifdef THREEDIMS
#define K_NORM (21.0 / (2.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define K_NORM (7.0 / M_PI) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define K_NORM (5.0 / 4.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL

#ifdef THREEDIMS
#define K_NORM (495.0 / (32.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define K_NORM (9.0 / M_PI) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define K_NORM (3.0 / 2.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL

#ifdef THREEDIMS
#define K_NORM (1365.0 / (64.0 * M_PI)) /*!< For 3D-normalized kernel */
#endif

#ifdef TWODIMS
#define K_NORM (78.0 / (7.0 * M_PI)) /*!< For 2D-normalized kernel */
#endif

#ifdef ONEDIMS
#define K_NORM (55.0 / 32.0) /*!< For 1D-normalized kernel */
#endif

#endif /* WENDLAND_C6_KERNEL */

static int int_compare(const void *a, const void *b);

/* sph loop kernel function -> u < 1 */
void kernel(double u, double hinv3, double hinv4, double *wk, double *dwk)
{
#ifdef CUBIC_SPLINE_KERNEL
#if defined(WENDLAND_C2_KERNEL) || defined(WENDLAND_C4_KERNEL) || defined(WENDLAND_C6_KERNEL)
#error "Only one SPH kernel can be used"
#endif
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
#endif

#ifdef WENDLAND_C2_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);

  
  *dwk = -12.0 * u * t2;
  
  *wk = t2 * t1 * (1.0 + u * 3.0);

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  
  *dwk = -20.0 * u * t2 * t1;
  
  *wk = t4 * (1.0 + u * 4.0);

#endif
#endif /* WENDLAND_C2_KERNEL */

#ifdef WENDLAND_C4_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = t1 * t1;
  double t4 = t2 * t2;
  double t5 = t4 * t1;

  
  *dwk = -14.0 * t4 * (4.0 * u + 1) * u;
  
  *wk = t5 * (1.0 + u * (5.0 + 8.0 * u));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t2 * t2 * t2;
  
  *dwk = -56.0 / 3.0 * u * t4 * t1 * (5.0 * u + 1);
  
  *wk = t6 * (1.0 + u * (6.0 + 35.0 / 3.0 * u));

#endif
#endif /* WENDLAND_C4_KERNEL */

#ifdef WENDLAND_C6_KERNEL /* Dehnen & Aly 2012 */
#ifdef ONEDIMS
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t6 = t4 * t2;
  double t7 = t4 * t2 * t1;
  
  *dwk = -6.0 * u * t6 * (3.0 + u * (18.0 + 35.0 * u));
  
  *wk = t7 * (1.0 + u * (7.0 + u * (19.0 + 21.0 * u)));

#else /* 2d or 3d */
  double t1 = (1.0 - u);
  double t2 = (t1 * t1);
  double t4 = t2 * t2;
  double t7 = t4 * t2 * t1;
  double t8 = t4 * t4;
  
  *dwk = -22.0 * u * (1.0 + u * (7.0 + 16.0 * u)) * t7;
  
  *wk = t8 * (1.0 + u * (8.0 + u * (25.0 + 32.0 * u)));

#endif
#endif /* WENDLAND_C6_KERNEL */
  
  *dwk *= K_NORM * hinv4;
  
  *wk *= K_NORM * hinv3;
}
/* THIS PART ADAPTED FROM GADGET4 */

#ifdef BLACKHOLES
#ifdef BONDI_ACCRETION 
void update_bh_accretion_rate(void)
{
  /* calculate bondi accretion rate */
  int i;
  double density, pressure, sound_speed, velocity_gas_norm;
  double denominator, denominator_inv, BondiRate, EddingtonRate;
  double accretion_rate, acc_rate_for_print;

  accretion_rate = acc_rate_for_print = 0;

  for(i = 0; i < NumBhs; i++)
    {
      /* get pressure */
      if(BhP[i].Density>0)
        {  
          density = BhP[i].Density;
          pressure = GAMMA_MINUS1 * density * BhP[i].InternalEnergyGas;

          /* get soundspeed */
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
  
      /* limit by Eddington accretion rate */
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

#ifdef STARS
/* update the SNII flag */
void update_SNII(void)
{
  int i;
  
  for(i=0; i<NumStars; i++)
    {
      if(SP[i].SNIIFlag == 1)
        {                         
          PPS(i).Mass = SP[i].SNIIRemnantMass;
          SP[i].SNIIFlag = 2;
          continue;
        }

      if(SP[i].SNIIFlag == 2)
        continue;

      if(All.Time > SP[i].SNIITime)
        {
          SP[i].SNIIFlag = 1; 
          
          //timebin_add_particle(&TimeBinsStar, i, -1, 0, 1);

          /* SNII feedback variables */
          double elements[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};

          struct CELibStructFeedbackInput Input = 
            {
              .Mass = (PPS(i).Mass * All.UnitMass_in_g / SOLAR_MASS),                 
              .Metallicity = 0.0004,          
              .MassConversionFactor = 1, 
              .Elements = elements,
            };

          struct CELibStructFeedbackOutput Output = 
            CELibGetFeedback(Input, CELibFeedbackType_SNII);

          SP[i].SNIIEnergyFeed = (Output.Energy / All.UnitEnergy_in_cgs);
          SP[i].SNIIMassFeed = (Output.EjectaMass * SOLAR_MASS / All.UnitMass_in_g);
          SP[i].SNIIRemnantMass = (Output.RemnantMass * SOLAR_MASS / All.UnitMass_in_g);
        }  
    } 
}
#endif

#ifdef BLACKHOLES
/* get timestep for bh based on smallest between ngbs */
integertime get_timestep_bh(int p)
{ 
  return BhP[p].NgbMinStep;
}

/* update bh-timestep at prior_mesh_construction */
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

/* call this function as the reconstruct_timebins() bh version */
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

/* call this function after updating the bh-timebin to the ngb condition */
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
#endif

#ifdef STARS
/* get timestep for star based on smallest between ngbs */
integertime get_timestep_star(int p)
{ 
  if(SP[p].SNIIFlag == 1)
    return 0;
  
  return (integertime)(1e-2 / All.Timebase_interval);
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

/* call this function as the reconstruct_timebins() star version */
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

/* call this function after updating the star-timebin to the ngb condition */
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
#endif

void perform_end_of_step_physics(void)
{
  int idx, i;
  double pj, p0;
  double kick_vector[3], momentum_kick[3];

  momentum_kick[0] = momentum_kick[1] = momentum_kick[2] = 0;

#ifdef BLACKHOLES
#ifdef BONDI_ACCRETION
  int j, bin;
  double dt;
  /* accrete mass, angular momentum onto the bh and drain ngb cells */
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
                  /* we're also losing thermal and kinetic energy & momentum */
                
                  /* update total energy */
                  SphP[j].Energy *= 0.1;

                  /* update momentum */
                  SphP[j].Momentum[0] *= 0.1;
                  SphP[j].Momentum[1] *= 0.1;
                  SphP[j].Momentum[2] *= 0.1;
                }
              else
                {
                  P[j].Mass -= SphP[j].MassDrain;
                
                  /* update total energy */
                  SphP[j].Energy *= (P[j].Mass)/(P[j].Mass + SphP[j].MassDrain);

                  /* update momentum */
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

    struct pv_update_data pvd;
    if(All.ComovingIntegrationOn)
      {
        pvd.atime    = All.Time;
        pvd.hubble_a = hubble_function(All.Time);
        pvd.a3inv    = 1 / (All.Time * All.Time * All.Time);
      }
    else
        pvd.atime = pvd.hubble_a = pvd.a3inv = 1.0;

    /* inject feedback to ngb cells */
    if(All.Time >= All.FeedbackTime)
    {   
      if(All.FeedbackFlag > 0)
        {
          for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
            {
              i = TimeBinsHydro.ActiveParticleList[idx];
              if(i < 0)
              continue;

#ifdef BLACKHOLES
              /* dump mass, momentum and energy injected by bh */
              if(SphP[i].ThermalFeed > 0 || SphP[i].KineticFeed > 0)
                {
                  /* add mass */
                  P[i].Mass += SphP[i].MassLoading;

                  /* add kinetic energy */
                  SphP[i].Energy += SphP[i].KineticFeed;

                  /* calculate momentum feed exactly so energy is conserved */
                  /*-> we need to do this here so that particle properties don't change between loading the buffer and emptying it*/
                  kick_vector[0] = SphP[i].BhKickVector[0];
                  kick_vector[1] = SphP[i].BhKickVector[1];
                  kick_vector[2] = SphP[i].BhKickVector[2];

                  p0 = sqrt(pow(SphP[i].Momentum[0], 2) + pow(SphP[i].Momentum[1], 2) + pow(SphP[i].Momentum[2], 2));
              
                  pj = sqrt(2 * P[i].Mass * (SphP[i].Energy - (P[i].Mass-SphP[i].MassLoading)*SphP[i].Utherm)) - p0;
    
                  momentum_kick[0] = kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
                  momentum_kick[1] = kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
                  momentum_kick[2] = kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2)); 
   
                  /* update total energy */
                  SphP[i].Energy += SphP[i].ThermalFeed;
                  All.EnergyExchange[1] += SphP[i].ThermalFeed + SphP[i].KineticFeed;
                  /* update momentum */
                  SphP[i].Momentum[0] += momentum_kick[0];
                  SphP[i].Momentum[1] += momentum_kick[1];
                  SphP[i].Momentum[2] += momentum_kick[2];
                  /* update velocities */
                  update_primitive_variables_single(P, SphP, i, &pvd);
                  /* update internal energy */
                  update_internal_energy(P, SphP, i, &pvd);
                  /* update pressure */
                  set_pressure_of_cell_internal(P, SphP, i);
                  /* set feed flags to zero */
                  SphP[i].ThermalFeed = SphP[i].KineticFeed = SphP[i].MassLoading = 0;
                  momentum_kick[0] = momentum_kick[1] = momentum_kick[2] = 0;
#ifdef PASSIVE_SCALARS                 
                  /* tracer field advected passively */
                  SphP[i].PScalars[0] = 1;
                  SphP[i].PConservedScalars[0] = P[i].Mass;
#endif
                }
#endif          
            }
#ifdef BURST_MODE
          All.FeedbackFlag = -1;
#endif 
        }    

      for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
        {
          i = TimeBinsHydro.ActiveParticleList[idx];
          if(i < 0)
          continue;

#ifdef STARS            
          /* dump mass, momentum and energy injected by stars */              
          if(SphP[i].MomentumFeed > 0 || SphP[i].ThermalEnergyFeed > 0 || SphP[i].KineticEnergyFeed > 0)
            {               
              /* add mass */
              //P[i].Mass += SphP[i].MassFeed;

              /* calculate kick */
              kick_vector[0] = SphP[i].MomentumKickVector[0];
              kick_vector[1] = SphP[i].MomentumKickVector[1];
              kick_vector[2] = SphP[i].MomentumKickVector[2];
              
              /******  momentum conserving wind *****/
              //pj = SphP[i].MomentumFeed;

              /* update momentum */
              //SphP[i].Momentum[0] += kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              //SphP[i].Momentum[1] += kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              //SphP[i].Momentum[2] += kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));  

              //All.EnergyExchange[3] += SphP[i].MomentumFeed;     
                 
              /* update velocities */
              //update_primitive_variables_single(P, SphP, i, &pvd);
              
              /***** energy conserving supernova *****/
              
              /* add kinetic energy */
              SphP[i].Energy += SphP[i].KineticEnergyFeed * All.cf_atime*All.cf_atime;

              /* calculate momentum feed exactly so energy is conserved */
              /*-> we need to do this here so that particle properties don't change between loading the buffer and emptying it*/
              p0 = sqrt(pow(SphP[i].Momentum[0], 2) + pow(SphP[i].Momentum[1], 2) + pow(SphP[i].Momentum[2], 2));

              pj = sqrt(2 * P[i].Mass * (SphP[i].Energy - (P[i].Mass/*-SphP[i].FMass*/)*SphP[i].Utherm /*- SphP[i].FMass * u*/)) - p0;

              momentum_kick[0] = kick_vector[0] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              momentum_kick[1] = kick_vector[1] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2));
              momentum_kick[2] = kick_vector[2] * pj / sqrt(pow(kick_vector[0], 2) + pow(kick_vector[1], 2) + pow(kick_vector[2], 2)); 

              /* update total energy */
              SphP[i].Energy += SphP[i].ThermalEnergyFeed * All.cf_atime*All.cf_atime;
              All.EnergyExchange[5] += SphP[i].ThermalEnergyFeed + SphP[i].KineticEnergyFeed; 
              /* update momentum */
              SphP[i].Momentum[0] += momentum_kick[0];
              SphP[i].Momentum[1] += momentum_kick[1];
              SphP[i].Momentum[2] += momentum_kick[2];
              /* update velocities */
              update_primitive_variables_single(P, SphP, i, &pvd);
              /* update internal energy */
              update_internal_energy(P, SphP, i, &pvd);
              /* update pressure */
              set_pressure_of_cell_internal(P, SphP, i);
              /* set feed flags to zero */
              SphP[i].MomentumFeed = SphP[i].ThermalEnergyFeed = SphP[i].KineticEnergyFeed = SphP[i].MassFeed = 0;
              momentum_kick[0] = momentum_kick[1] = momentum_kick[2] = 0;
            }
#endif
        }
       
    }
  MPI_Allreduce(&All.EnergyExchange, &All.EnergyExchangeTot, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // synchronize all tasks
  
#ifdef BLACKHOLES
  mpi_printf("BLACK_HOLES: Energy given by BH = %e, Energy taken up by gas particles = %e \n", 
    All.EnergyExchangeTot[0] * All.UnitEnergy_in_cgs, All.EnergyExchangeTot[1] * All.UnitEnergy_in_cgs);
#endif /* ifdef BLACKHOLES */
  mpi_printf("STARS: Momentum given by StarParts = %e, Momentum taken up by gas particles = %e \n", 
    All.EnergyExchangeTot[2] * All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s, 
    All.EnergyExchangeTot[3] * All.UnitMass_in_g * All.UnitVelocity_in_cm_per_s);
  mpi_printf("STARS: Energy given by StarParts = %e, Energy taken up by gas particles = %e \n", 
    All.EnergyExchangeTot[4] * All.UnitEnergy_in_cgs, All.EnergyExchangeTot[5] * All.UnitEnergy_in_cgs);

#ifdef BURST_MODE
  if(All.EnergyExchangeTot[0] - All.EnergyExchangeTot[1] > 10)  
    All.FeedbackFlag = 1;
#endif  
}

static int int_compare(const void *a, const void *b)
{
  if(*((int *)a) < *((int *)b))
    return -1;

  if(*((int *)a) > *((int *)b))
    return +1;

  return 0;
}
#endif
