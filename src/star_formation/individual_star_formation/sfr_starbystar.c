#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../gravity/forcetree.h"

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
  double dens, temp;

  /* note: assuming FULL ionization */ //need grackle fields to do properly
  double u_to_temp_fac =
  (4 / (8 - 5 * (1 - HYDROGEN_MASSFRAC))) * PROTONMASS / BOLTZMANN * GAMMA_MINUS1 * All.UnitEnergy_in_cgs / All.UnitMass_in_g;

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

      dens = SphP[i].Density;
      temp = SphP[i].Utherm * u_to_temp_fac;

      double numberdens_threshold = All.NumberDensThreshold * PROTONMASS / All.UnitDensity_in_cgs;

      /* check whether conditions for star formation are fulfilled.
       * f=1  normal cooling
       * f=0  star formation
       */

      flag = 1; /* default is normal cooling */

      /* enable star formation if gas is above SF density threshold */
      if(dens * All.cf_a3inv >= numberdens_threshold && temp < All.TemperatureThreshold)
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