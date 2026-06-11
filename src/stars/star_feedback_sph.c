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
  MyDouble NgbsMass;
  MyDouble NgbsVolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeSN_yr;
  MyDouble PhysicalAge_yr;
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

  MyFloat Hsml;
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
  for(int j = 0; j < 3; j++)
    {  
      in->Pos[j] = PPS(i).Pos[j];
      in->Vel[j] = PPS(i).Vel[j];
    }
  in->NgbsMass = SP[i].NgbsMass;
  in->NgbsVolume = SP[i].NgbsVolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  in->TimeSN_yr         = SP[i].TimeSN_yr;
  in->PhysicalAge_yr      = SP[i].PhysicalAge_yr;
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

  in->Hsml = SP[i].Hsml;
  in->Firstnode = firstnode;
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
  TIMER_START(CPU_STARS_FEEDBACK);
  
  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsStar.NActiveParticles, kernel_local, kernel_imported);

  TIMER_STOP(CPU_STARS_FEEDBACK);
}

static int star_feedback_evaluate(int target, int mode, int threadid)
{
  int j, n, numnodes, *firstnode; 
  double h, h2, dx, dy, dz, r, r2, wk; 
  MyDouble *pos, *vel, ngbsmass, ngbsvolume, factor;
  
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
  
  ngbsmass = target_data->NgbsMass;
  ngbsvolume = target_data->NgbsVolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble time_sn_yr = target_data->TimeSN_yr;
  MyDouble physical_age_yr = target_data->PhysicalAge_yr;
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

      if(P[j].Type != 0 || P[j].Mass == 0 || P[j].ID == 0)
        continue;

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

          factor = SphP[j].Volume / ngbsvolume;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
          if(time_sn_yr < MAX_REAL_NUMBER)
            {
              double E_inject_code = 1e51 / 
              (All.cf_UnitMass_in_g * All.cf_UnitVelocity_in_cm_per_s * All.cf_UnitVelocity_in_cm_per_s);

              double unew = SphP[j].Utherm + E_inject_code * factor / P[j].Mass;

              double t_frac = physical_age_yr / time_sn_yr;
              t_frac = fmin(fmax(t_frac, 0.0), 1.0);

              double Csn = SphP[j].Csnd + (sqrt(GAMMA * GAMMA_MINUS1 * unew) - SphP[j].Csnd) * t_frac;
          
              if(Csn > SphP[j].Csn)
                SphP[j].Csn = Csn;
            }
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
              SphP[j].StarMomentumFeed[0] += (windmomentum * dx/r + massloss * vel[0] / All.cf_atime) * factor;
              SphP[j].StarMomentumFeed[1] += (windmomentum * dy/r + massloss * vel[1] / All.cf_atime) * factor;
              SphP[j].StarMomentumFeed[2] += (windmomentum * dz/r + massloss * vel[2] / All.cf_atime) * factor;
              
              double vsq = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
              double vdotp = vel[0] * (dx / r) + vel[1] * (dy / r) + vel[2] * (dz / r);
              SphP[j].StarEnergyFeed += (windmomentum * windmomentum / (2.0 * massloss) // wind KE
                                     + 0.5 * massloss * vsq / All.cf_atime / All.cf_atime       // bulk KE
                                     + windmomentum * vdotp / All.cf_atime) * factor;       // cross
              All.StarFeedbackLocal[2] += (windmomentum * windmomentum / (2.0 * massloss) 
                                       + 0.5 * massloss * vsq / All.cf_atime / All.cf_atime                             
                                       + windmomentum * vdotp / All.cf_atime) * factor;
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

              SphP[j].StarMomentumFeed[0] += (pSN * dx/r + SNmassloss * vel[0] / All.cf_atime) * factor;
              SphP[j].StarMomentumFeed[1] += (pSN * dy/r + SNmassloss * vel[1] / All.cf_atime) * factor;
              SphP[j].StarMomentumFeed[2] += (pSN * dz/r + SNmassloss * vel[2] / All.cf_atime) * factor;

              SphP[j].StarEnergyFeed += (SNenergyinject                              // thermal/KE in star frame
                                     + 0.5 * SNmassloss * vsq / All.cf_atime / All.cf_atime  // bulk KE
                                     + pSN * vdotp / All.cf_atime) * factor;             // cross
              All.StarFeedbackLocal[2] += (SNenergyinject       
                                       + 0.5 * SNmassloss * vsq / All.cf_atime / All.cf_atime  
                                       + pSN * vdotp / All.cf_atime) * factor;     
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

              // Momentum in sim frame = bulk motion of ejected mass only 
              SphP[j].StarMomentumFeed[0] += SNmassloss * vel[0] / All.cf_atime * factor;
              SphP[j].StarMomentumFeed[1] += SNmassloss * vel[1] / All.cf_atime * factor;
              SphP[j].StarMomentumFeed[2] += SNmassloss * vel[2] / All.cf_atime * factor;

              // Energy = thermal (in star frame) + bulk KE of ejected mass (no cross term)
              SphP[j].StarEnergyFeed += (SNenergyinject
                                     + 0.5 * SNmassloss * vsq / All.cf_atime / All.cf_atime) * factor;
              All.StarFeedbackLocal[2] += (SNenergyinject
                                       + 0.5 * SNmassloss * vsq / All.cf_atime / All.cf_atime) * factor;
            }
#endif
        }
    }

  return 0;
}