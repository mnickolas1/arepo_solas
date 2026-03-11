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
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyFloat Hsml;
  MyDouble NgbMass;
  MyDouble NgbVolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeSN;
  MyDouble Birthtime;
#endif  

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
  in->Pos[0]         = PPS(i).Pos[0];
  in->Pos[1]         = PPS(i).Pos[1];
  in->Pos[2]         = PPS(i).Pos[2];
  in->Vel[0]         = PPS(i).Vel[0];
  in->Vel[1]         = PPS(i).Vel[1];
  in->Vel[2]         = PPS(i).Vel[2];
  in->Hsml           = SP[i].Hsml;
  in->NgbMass        = SP[i].NgbMass;
  in->NgbVolume      = SP[i].NgbVolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  in->TimeSN         = SP[i].TimeSN;
  in->Birthtime      = SP[i].Birthtime;
#endif  

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
  int idx, i, j;

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
      
      if(SP[i].Active == 1)    
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
  double h, h2, dx, dy, dz, r, r2, wk; 
  MyDouble *pos, *vel, ngbmass, ngbvolume, factor;
  
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
  
  pos = target_data->Pos;
  vel = target_data->Vel;
  h = target_data->Hsml;
  h2   = h * h;
  
  ngbmass = target_data->NgbMass;
  ngbvolume = target_data->NgbVolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble timesn = target_data->TimeSN;
  MyDouble birthtime = target_data->Birthtime;
#endif

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

          factor = SphP[j].Volume / ngbvolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
          double unew = (SphP[j].Utherm * P[j].Mass + (1e51 / All.UnitEnergy_in_cgs) * factor) / (P[j].Mass /*+dm_inject*/);

          double t_frac = (All.Time - birthtime) / (timesn - birthtime);
          t_frac = fmin(fmax(t_frac, 0.0), 1.0);

          if(timesn < MAX_REAL_NUMBER)
            SphP[j].Csn = SphP[j].Csnd + (sqrt(GAMMA * GAMMA_MINUS1 * unew) - SphP[j].Csnd) * t_frac;
#endif
              
#ifdef WINDS
          if(massloss > 0)
            {  
              SphP[j].StarMassFeed += massloss * factor;
              All.StarFeedbackLocal[0] += massloss * factor;
#ifdef METALS
              SphP[j].StarMetalsFeed += metalsloss * factor;
              All.StarFeedbackLocal[1] += metalsloss * factor;
#endif
              SphP[j].StarMomentumFeed[0] += (windmomentum * dx/r + massloss * vel[0]) * factor;
              SphP[j].StarMomentumFeed[1] += (windmomentum * dy/r + massloss * vel[1]) * factor;
              SphP[j].StarMomentumFeed[2] += (windmomentum * dz/r + massloss * vel[2]) * factor;
              All.StarFeedbackLocal[2] += windmomentum * factor; //need to update 
              
              double vsq = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
              double vdotp = vel[0] * (dx / r) + vel[1] * (dy / r) + vel[2] * (dz / r);
              SphP[j].StarEnergyFeed += (windmomentum * windmomentum / (2.0 * massloss) // wind KE
                                     + 0.5 * massloss * vsq                             // bulk KE
                                     + windmomentum * vdotp) * factor;                     // cross
              All.StarFeedbackLocal[3] += (windmomentum * windmomentum / (2.0 * massloss) 
                                       + 0.5 * massloss * vsq                             
                                       + windmomentum * vdotp) * factor;
            }
#endif     

#ifdef SUPERNOVAE
          if(SNmassloss > 0)
            {
              SphP[j].StarMassFeed += SNmassloss * factor;
              All.StarFeedbackLocal[0] += SNmassloss * factor;
#ifdef METALS
              SphP[j].StarMetalsFeed += SNmetalsloss * factor;
              All.StarFeedbackLocal[1] += SNmetalsloss * factor;
#endif

              double pSN = sqrt(2.0 * SNenergyinject * SNmassloss);  // |p| in star frame
              double vsq = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
              double vdotp = vel[0]*(dx/r) + vel[1]*(dy/r) + vel[2]*(dz/r);

              SphP[j].StarMomentumFeed[0] += (pSN * dx/r + SNmassloss * vel[0]) * factor;
              SphP[j].StarMomentumFeed[1] += (pSN * dy/r + SNmassloss * vel[1]) * factor;
              SphP[j].StarMomentumFeed[2] += (pSN * dz/r + SNmassloss * vel[2]) * factor;
              All.StarFeedbackLocal[2] += sqrt(2 * SNenergyinject * SNmassloss) * factor; // need to update

              SphP[j].StarEnergyFeed += (SNenergyinject        // thermal/KE in star frame
                                     + 0.5 * SNmassloss * vsq  // bulk KE
                                     + pSN * vdotp) * factor;     // cross
              All.StarFeedbackLocal[3] += (SNenergyinject       
                                       + 0.5 * SNmassloss * vsq  
                                       + pSN * vdotp) * factor;     
            }
#endif

#ifdef SUPERNOVAE_THERMAL
          if(SNmassloss > 0)
            {
              SphP[j].StarMassFeed += SNmassloss * factor;
              All.StarFeedbackLocal[0] += SNmassloss * factor;
#ifdef METALS
              SphP[j].StarMetalsFeed += SNmetalsloss * factor;
              All.StarFeedbackLocal[1] += SNmetalsloss * factor;
#endif
              double vsq = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];

              // Momentum in sim frame = bulk motion of ejected mass only (no r̂ term)
              SphP[j].StarMomentumFeed[0] += SNmassloss * vel[0] * factor;
              SphP[j].StarMomentumFeed[1] += SNmassloss * vel[1] * factor;
              SphP[j].StarMomentumFeed[2] += SNmassloss * vel[2] * factor;
              All.StarFeedbackLocal[2] += 0.0; // zero momentum injected in star frame

              // Energy = thermal (in star frame) + bulk KE of ejected mass (no cross term)
              SphP[j].StarEnergyFeed += (SNenergyinject
                                     + 0.5 * SNmassloss * vsq) * factor;
              All.StarFeedbackLocal[3] += (SNenergyinject
                                       + 0.5 * SNmassloss * vsq) * factor;
            }
#endif
        }
    }

  return 0;
}