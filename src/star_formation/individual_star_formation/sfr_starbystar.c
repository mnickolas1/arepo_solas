#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

static double compute_mu(int i)
{
    /* Grab mass fractions from Grackle */
    double XH  = SphP[i].grHI + SphP[i].grHII; // atomic hydrogen fraction
    double XH2 = SphP[i].grH2I; // molecular hydrogen fraction 
    double XHe = SphP[i].grHeI + SphP[i].grHeII + SphP[i].grHeIII; // helium fraction
    double Xe  = SphP[i].Ne; // electron fraction

    /* Optional: include metals if desired (usually negligible) */
    double Z = SphP[i].Metallicity;
    double A_Z = 16.0 * PROTONMASS; // metals, approx oxygen

    /* Compute mean molecular weight: g per particle */
    /* 1 H atom = m_H, 1 He atom = 4*m_H, electrons = m_H (counted as particles for n) */
    return 1.0 / (XH + XH2/2.0 + XHe/4.0 + Xe + Z/16.0); // -> Dimensionless mu
}

/*! \brief Main driver for star formation and gas cooling.
 *
 *
 *  \return void
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i, flag;
  double dt, du, unew;
  double number_dens, temp;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      if(P[i].Mass == 0 && P[i].ID == 0)
        continue; /* skip cells that have been swallowed or eliminated */

      dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      dt *= All.cf_atime  / All.cf_time_hubble_a;

      /* apply the temperature floor */
      unew = dmax(All.MinEgySpec, SphP[i].Utherm);

      if(unew < 0)
        terminate("Invalid Temperature: Task=%d i=%d unew=%g\n", ThisTask, i, unew);

      du = unew - SphP[i].Utherm;
      SphP[i].Utherm += du;
      SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

      cool_cell(i);
      
      //double mu = compute_mu(i); 
      double mu = 2.33; // molecular H

      number_dens = (SphP[i].Density * All.UnitDensity_in_cgs) / mu / PROTONMASS;
      number_dens /= All.cf_a3inv;
      
      double u_to_temp_fac = mu * PROTONMASS / BOLTZMANN * GAMMA_MINUS1;

      temp = (SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g) * u_to_temp_fac;
     
      /* check whether conditions for star formation are fulfilled.
       * f=1  normal cooling
       * f=0  star formation
       */

      flag = 1; /* default is normal cooling */

      /* enable star formation if gas is above SF density threshold */
      if(number_dens >= All.NumberDensThreshold && temp < All.TemperatureThreshold)
        if(All.Time > 0)  
          flag = 0;

      /* tracer particles don't form stars */    
      if(P[i].Mass == 0) 
        flag = 1;
      
      /* inactive star formation */
      if(flag == 1)
          SphP[i].Sfr = 0;

      /* active star formation */
      if(flag == 0)
        {
          SphP[i].Sfr = 1;
        }

    }
  TIMER_STOP(CPU_COOLINGSFR);
}