#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

static int star_feedback_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  int Bin;  
  MyDouble Pos[3];
  MyFloat Hsml;
  MyDouble NgbMass;
  MyDouble NgbVolume;

#ifdef WINDS
  MyDouble MassLoss;
#ifdef METALS
  MyDouble MetalsLoss;
#endif
  MyDouble WindMomentum;
#endif

#ifdef SUPERNOVAE
  MyDouble SN_MassLoss;
#ifdef METALS
  MyDouble SN_MetalsLoss;
#endif
  MyDouble SN_EnergyInject;
#endif

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
  in->NgbMass        = SP[i].NgbMass;
  in->NgbVolume      = SP[i].NgbVolume;

#ifdef WINDS
  in->MassLoss       = SP[i].MassLoss;
#ifdef METALS
  in->MetalsLoss     = SP[i].MetalsLoss;
#endif
  in->WindMomentum   = SP[i].WindMomentum;
#endif

#ifdef SUPERNOVAE
  in->SN_MassLoss    = SP[i].SN_MassLoss;
#ifdef METALS
  in->SN_MetalsLoss  = SP[i].SN_MetalsLoss;
#endif
  in->SN_EnergyInject = SP[i].SN_EnergyInject;
#endif

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
  int i, idx, j;

  int threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      //i = NextParticle++;

      //if(i >= NumStars)
      //  break;
        
      idx = NextParticle++;

      if(idx >= TimeBinsStar.NActiveParticles)
        break;

      i = TimeBinsStar.ActiveParticleList[idx];

      star_feedback_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

  int threadid = get_thread_num();

  while(1)
    {
      i = cnt++;

      if(i >= Nimport)
        break;

      star_feedback_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

void star_feedback(void)
{
  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsStar.NActiveParticles, kernel_local, kernel_imported);
}

static int star_feedback_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode; 
  int bin;
  double h, h2, hinv, hinv3, hinv4; 
  double dx, dy, dz, r, r2, u, wk, dwk;
  MyDouble *pos, ngbmass, ngbvolume;
  
  data_in local, *target_data;

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
  
  bin       = target_data->Bin;
  pos       = target_data->Pos;
  h         = target_data->Hsml;
  
  ngbmass   = target_data->NgbMass;
  ngbvolume = target_data->NgbVolume;

#ifdef WINDS
  MyDouble massloss = target_data->MassLoss;
#ifdef METALS
  MyDouble metalsloss = target_data->MetalsLoss;
#endif
  MyDouble windmomentum = target_data->WindMomentum;
#endif

#ifdef SUPERNOVAE
  MyDouble SNmassloss = target_data->SN_MassLoss;
#ifdef METALS
  MyDouble SNmetalsloss = target_data->SN_MetalsLoss;
#endif
  MyDouble SNenergyinject = target_data->SN_EnergyInject;
#endif  

  h2   = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;
 
/* star timestep */
  double dt = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;
  dt *= All.cf_atime / All.cf_time_hubble_a;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);
  for(n = 0; n < nfound; n++)
    {
      j = Thread[threadid].Ngblist[n];

      dx = P[j].Pos[0] - pos[0];
      dy = P[j].Pos[1] - pos[1]; 
      dz = P[j].Pos[2] - pos[2]; 

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

          star_kernel(u, hinv3, hinv4, &wk, &dwk);
              
#ifdef WINDS
          if(massloss > 0)
            {  
              SphP[j].StarMassFeed += massloss * SphP[j].Volume / ngbvolume;
              All.StarFeedbackLocal[0] += massloss * SphP[j].Volume / ngbvolume;
#ifdef METALS
              SphP[j].StarMetalsFeed += metalsloss * SphP[j].Volume / ngbvolume;
              All.StarFeedbackLocal[1] += metalsloss * SphP[j].Volume / ngbvolume;
#endif
          
              SphP[j].StarMomentumFeed[0] += windmomentum * dx/r * SphP[j].Volume / ngbvolume;
              SphP[j].StarMomentumFeed[1] += windmomentum * dy/r * SphP[j].Volume / ngbvolume;
              SphP[j].StarMomentumFeed[2] += windmomentum * dz/r * SphP[j].Volume / ngbvolume;
              All.StarFeedbackLocal[2] += windmomentum * SphP[j].Volume / ngbvolume;
          
              SphP[j].StarEnergyFeed += windmomentum *  windmomentum / (2 * massloss) * SphP[j].Volume / ngbvolume;
              All.StarFeedbackLocal[3] += windmomentum *  windmomentum / (2 * massloss) * SphP[j].Volume / ngbvolume;
            }
#endif     

#ifdef SUPERNOVAE
          if(SNmassloss > 0)
            {
              SphP[j].StarMassFeed += SNmassloss * SphP[j].Volume / ngbvolume;
              All.StarFeedbackLocal[0] += SNmassloss * SphP[j].Volume / ngbvolume;
#ifdef METALS
              SphP[j].StarMetalsFeed += SNmetalsloss * SphP[j].Volume / ngbvolume;
              All.StarFeedbackLocal[1] += SNmetalsloss * SphP[j].Volume / ngbvolume;
#endif
              SphP[j].StarMomentumFeed[0] += sqrt(2 * SNenergyinject * SNmassloss) * dx/r * SphP[j].Volume / ngbvolume;
              SphP[j].StarMomentumFeed[1] += sqrt(2 * SNenergyinject * SNmassloss) * dy/r * SphP[j].Volume / ngbvolume;
              SphP[j].StarMomentumFeed[2] += sqrt(2 * SNenergyinject * SNmassloss) * dz/r * SphP[j].Volume / ngbvolume;   
              All.StarFeedbackLocal[2] += sqrt(2 * SNenergyinject * SNmassloss) * SphP[j].Volume / ngbvolume;

              SphP[j].StarEnergyFeed += SNenergyinject * SphP[j].Volume / ngbvolume;
              All.StarFeedbackLocal[3] += SNenergyinject * SphP[j].Volume / ngbvolume;
            }
#endif
        }
    }

  return 0;
}