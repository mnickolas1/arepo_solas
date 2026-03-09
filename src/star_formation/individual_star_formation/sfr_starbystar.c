#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../main/allvars.h"
#include "../../main/proto.h"

/*! \brief Main driver for star formation and gas cooling.
 *
 *
 *  \return void
 */
void cooling_and_starformation(void)
{
  TIMER_START(CPU_COOLINGSFR);

  int idx, i;
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

      cool_cell(i); //do we need another temp floor?
      
      //double mu = compute_mu(i); 
      double mu = 2.33; // molecular H

      number_dens = (SphP[i].Density * All.UnitDensity_in_cgs) / mu / PROTONMASS;
      //number_dens /= All.cf_a3inv; //need to figure out a factors
      
      double u_to_temp_fac = mu * PROTONMASS / BOLTZMANN * GAMMA_MINUS1;
      temp = (SphP[i].Utherm * All.UnitEnergy_in_cgs / All.UnitMass_in_g) * u_to_temp_fac;
      
      /* default is just cooling */
      SphP[i].Sfr = 0; 

      /* star formation if gas is dense and cold */
      if(number_dens >= All.NumberDensThreshold && temp < All.TemperatureThreshold)
        if(All.Time > 0)  
          SphP[i].Sfr = 1;

      /* tracer particles don't form stars */    
      if(P[i].Mass == 0) 
        SphP[i].Sfr = 0;
    }
  TIMER_STOP(CPU_COOLINGSFR);
}