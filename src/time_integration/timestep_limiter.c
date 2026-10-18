/*! \file timestep_limiter.c
 *  \brief Saitoh & Makino style neighbour time-bin limiter.
 *
 *  Determines, for every active hydro cell, the minimum hydro timebin among its
 *  Voronoi connection neighbours, so that find_timesteps_without_gravity() can cap
 *  the cell's step
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"
#include "../mesh/voronoi/voronoi.h"

/*! \brief Minimum connection-neighbour hydro time bin of every active cell
 *
 *  \return void
 */
void timestep_limiter_min_ngb_bin(int *min_ngb_bin)
{
  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
     /* No constraint */
      min_ngb_bin[idx] = TIMEBINS; 

      int i = TimeBinsHydro.ActiveParticleList[idx];

      if(i < 0)
        continue;

      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("Strange connectivity q=%d MaxNvc=%d", q, MaxNvc);

          int dp = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;

          if(particle >= 0)
            {
              int bin;

              if(DC[q].task != ThisTask)
                bin = PrimExch[particle].TimeBinHydro;
              else if(particle >= NumGas) 
                bin = P[particle - NumGas].TimeBinHydro;
              else
                bin = P[particle].TimeBinHydro;

              if(bin > 0 && bin < min_ngb_bin[idx])
                min_ngb_bin[idx] = bin;
            }

          if(q == SphP[i].last_connection)
          break;

          q = DC[q].next;
        }
    }
}