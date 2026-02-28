#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../extern/chealpix.h"
#include "../stars/star_radiation.h"

double HealpixDirs[MAX_RAYS][3];
int NRays; // 12 * NSIDE^2

double RayColumnDensity[MAX_RAYS];
double RayIntensity[MAX_RAYS];

void init_healpix_rays(int nside)
{
    NRays = 12 * nside * nside;
    for(int ipix=0; ipix<NRays; ipix++)
    {
        pix2vec_nest(nside, ipix, HealpixDirs[ipix]);
    }
}

void trace_ray(int star_index, double *dir, double *column_density, double *intensity)
{
  int i, idx;
   
  for(idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      i = TimeBinsStar.ActiveParticleList[idx];
    
    // Loop over rays
      for(int iray=0; iray<NRays; iray++)
        {
          double *dir = HealpixDirs[iray];

        // Initialize per-ray quantity
        RayColumnDensity[iray] = 0;
        RayIntensity[iray] = SP[i].RAD_Ionizing;

        // Walk the tree for this ray
        //raytrace_tree_walk(pos, dir, &RayColumnDensity[iray], &RayIntensity[iray]);
        }

    // Store results for this star, e.g. total absorbed energy or column
    }
}