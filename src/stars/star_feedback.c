#include <gsl/gsl_math.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"


static int star_ngb_feedback_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  int Bin;  
  MyDouble Pos[3];
  MyFloat Hsml;
  MyDouble StarMass;
  MyDouble StarDensity;
  MyDouble NgbMass;
  MyDouble NgbVolume;
  int SNIIFlag;
  MyDouble SNIIEnergyFeed;
  MyDouble SNIIMassFeed;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  in->Bin            = SP[i].TimeBinStar;
  in->Pos[0]         = PPS(i).Pos[0];
  in->Pos[1]         = PPS(i).Pos[1];
  in->Pos[2]         = PPS(i).Pos[2];
  in->Hsml           = SP[i].Hsml;
  in->StarMass       = PPS(i).Mass;
  in->StarDensity    = SP[i].Density;
  in->NgbMass        = SP[i].NgbMass;
  in->NgbVolume      = SP[i].NgbVolume;
  in->SNIIFlag       = SP[i].SNIIFlag;
  in->SNIIEnergyFeed = SP[i].SNIIEnergyFeed;
  in->SNIIMassFeed   = SP[i].SNIIMassFeed;
  in->Firstnode      = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.*/

typedef struct
{
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void*/
 
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
    }
  else /* combine */
    {
    }
}

#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int idx;

  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        idx = NextParticle++;

        if(idx >= TimeBinsStar.NActiveParticles)
          break;

        int i = TimeBinsStar.ActiveParticleList[idx];
        if(i < 0)
          continue;
        star_ngb_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        star_ngb_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

void star_ngb_feedback(void)
{
  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsStar.NActiveParticles, kernel_local, kernel_imported);
}

static int star_ngb_feedback_evaluate(int target, int mode, int threadid)
{
  int j, n, bin, snIIflag; 
  int numnodes, *firstnode;
  double h, h2, hinv, hinv3, hinv4; 
  double dx, dy, dz, r, r2, u, wk, dwk, dt;
  MyDouble *pos, star_mass, star_density, ngbmass, ngbvolume, snIImassfeed, snIIenergyfeed;

  data_in local, *target_data;
  //data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }
  
  bin            = target_data->Bin;
  pos            = target_data->Pos;
  h              = target_data->Hsml;
  star_mass      = target_data->StarMass;
  star_density   = target_data->StarDensity;
  ngbmass        = target_data->NgbMass;
  ngbvolume      = target_data->NgbVolume;
  snIIflag       = target_data->SNIIFlag;
  snIIenergyfeed = target_data->SNIIEnergyFeed;
  snIImassfeed   = target_data->SNIIMassFeed;

  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;
 
/* star timestep */
  //dt = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;
  dt  = All.TimeStep;
  dt *= All.cf_atime / All.cf_time_hubble_a;

  /* stellar wind */    
  double EddingtonLuminosity = 4. * M_PI * GRAVITY * (star_mass * All.UnitMass_in_g) * PROTONMASS * CLIGHT / THOMPSON;
  EddingtonLuminosity *=  (All.UnitTime_in_s / (All.UnitMass_in_g*pow(All.UnitVelocity_in_cm_per_s,2)));
  double energyfeed = EddingtonLuminosity * dt;
  /* supernova */    
  //if(snIIflag > 0)
  //  energyfeed = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);
  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = pos[0] - P[j].Pos[0];
      dy = pos[1] - P[j].Pos[1];
      dz = pos[2] - P[j].Pos[2];

#ifndef REFLECTIVE_X
      if(dx > boxHalf_X)
        dx -= boxSize_X;
      if(dx < -boxHalf_X)
        dx += boxSize_X;
#endif /* #ifndef REFLECTIVE_X */

#ifndef REFLECTIVE_Y
      if(dy > boxHalf_Y)
        dy -= boxSize_Y;
      if(dy < -boxHalf_Y)
        dy += boxSize_Y;
#endif /* #ifndef REFLECTIVE_Y */

#ifndef REFLECTIVE_Z
      if(dz > boxHalf_Z)
        dz -= boxSize_Z;
      if(dz < -boxHalf_Z)
        dz += boxSize_Z;
#endif /* #ifndef REFLECTIVE_Z */
      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2)
        {
          r = sqrt(r2);

          u = r * hinv;

          kernel(u, hinv3, hinv4, &wk, &dwk);

          /* set radial momentum kick for wind */
          /* uncomment for kernel */ 
          //SphP[j].MomentumFeed  += All.Lambda * energyfeed / (CLIGHT / All.UnitVelocity_in_cm_per_s) * P[j].Mass / star_density * wk;
          //All.EnergyExchange[2] += All.Lambda * energyfeed / (CLIGHT / All.UnitVelocity_in_cm_per_s) * P[j].Mass / star_density * wk;

          SphP[j].MomentumFeed  += All.Lambda * energyfeed / (CLIGHT / All.UnitVelocity_in_cm_per_s) * SphP[j].Volume / ngbvolume;
          All.EnergyExchange[2] += All.Lambda * energyfeed / (CLIGHT / All.UnitVelocity_in_cm_per_s) * SphP[j].Volume / ngbvolume;
              
          SphP[j].MomentumKickVector[0] = -dx;
          SphP[j].MomentumKickVector[1] = -dy;
          SphP[j].MomentumKickVector[2] = -dz;

          /* do supernova */
          if (snIIflag == 1)
            {              
              SphP[j].ThermalEnergyFeed += All.Fsn*All.Ftherm*snIIenergyfeed * SphP[j].Volume / ngbvolume;
              SphP[j].KineticEnergyFeed += All.Fsn*(1-All.Ftherm)*snIIenergyfeed * SphP[j].Volume / ngbvolume;
              
              All.EnergyExchange[4] += All.Fsn*snIIenergyfeed * SphP[j].Volume / ngbvolume;
              
              //SphP[j].MassFeed += snIImassfeed * SphP[j].Volume / ngbvolume;
            }
        }
    }

   /* Now collect the result at the right place */
  /*if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;*/

  return 0;
}

