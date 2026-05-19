#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

/* Pass counter: 1 = find host cell, 2 = gather feedback properties */
static int pass;

static int star_density_evaluate(int target, int mode, int threadid);
static int star_density_isactive(int n);

static MyFloat *StarNgbs;
struct HostCell 
{
  MyDouble Distance;
  int Index;
  int Task;
} *StarHostCell;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  
  /* Pass 2 inputs */
  struct HostCell HostCell;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  MyDouble TimeSN_yr;
  MyDouble PhysicalAge_yr;
#endif  

#if defined(WINDS) || defined(SUPERNOVAE)
  Mechanical_Feedback WindsAndSN;
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
    in->Pos[j] = PPS(i).Pos[j];

  /* Pass 2 inputs */  
  in->HostCell.Distance = StarHostCell[i].Distance;
  in->HostCell.Index = StarHostCell[i].Index;
  in->HostCell.Task = StarHostCell[i].Task;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
  in->TimeSN_yr = SP[i].TimeSN_yr;
  in->PhysicalAge_yr = SP[i].PhysicalAge_yr;
#endif  

#if defined(WINDS) || defined(SUPERNOVAE)
  in->WindsAndSN = SP[i].WindsAndSN;
#endif

  in->Hsml = SP[i].Hsml;
  in->Firstnode = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{ 
  /* Pass 1 outputs */
  MyDouble Ngbs;
  struct HostCell HostCell;

  /* Pass 2 outputs */
  int HostHydroBin;
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
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      /* Pass 1 outputs */
      if(pass == 1)
        {
          StarNgbs[i] = out->Ngbs;
          StarHostCell[i].Distance = out->HostCell.Distance;
          StarHostCell[i].Index = out->HostCell.Index;
          StarHostCell[i].Task = out->HostCell.Task;
        }

      /* Pass 2 outputs */
      if(pass == 2)
        {
          SP[i].HostHydroBin = out->HostHydroBin;
        }
    }
  else /* combine */
    {
      /* Pass 1 outputs */
      if(pass == 1)
        {
          StarNgbs[i] += out->Ngbs;
          if(out->HostCell.Distance < StarHostCell[i].Distance)
            {
              StarHostCell[i].Distance = out->HostCell.Distance;
              StarHostCell[i].Index = out->HostCell.Index;
              StarHostCell[i].Task = out->HostCell.Task;
            }
        }

      /* Pass 2 outputs */
      if(pass == 2)
        {
          if(!SP[i].HostHydroBin)
            SP[i].HostHydroBin = out->HostHydroBin;
        }
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
      
      double star_mass = PPS(i).Mass * All.cf_UnitMass_in_Msun;
      
      if(star_mass > 2)
        if(star_density_isactive(i))
          star_density_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
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

      star_density_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

/*! \brief Main function of SPH density calculation.
 *
 *  This function computes the local density for each active SPH particle and
 *  the number of weighted neighbors in the current smoothing radius. If a
 *  particle with its smoothing region is fully inside the local domain, it is
 *  not exported to the other processors. The function also detects particles
 *  that have a number of neighbors outside the allowed tolerance range. For
 *  these particles, the smoothing length is adjusted accordingly, and the
 *  computation is called again.
 *
 *  \return void
 */
void star_density(void)
{
  int idx, i, npleft, iter = 0;
  long long ntot;
  double t0, t1;

  CPU_Step[CPU_MISC] += measure_time();

  pass = 0;

  StarNgbs = (MyFloat *)mymalloc("StarNgbs", NumStars * sizeof(MyFloat));
  StarHostCell = (struct HostCell *)mymalloc("StarHostCell", NumStars * sizeof(struct HostCell));

  for(i = 0; i < NumStars; i++)
    {
      SP[i].DensityFlag = 1;
      
      if(SP[i].Hsml <= 0)
        SP[i].Hsml = All.SofteningTable[PPS(i).SofteningType]; // is this comoving?
    }

  generic_set_MaxNexport();
  
  /* 
  Pass 1 — Expand Hsml until we enclose at least one gas cell,
           then record the closest cell as the host.
  */

  pass++;

  do 
    {
      t0 = second();

      generic_comm_pattern(TimeBinsStar.NActiveParticles, kernel_local, kernel_imported);

      for(idx = 0, npleft = 0; idx < TimeBinsStar.NActiveParticles; idx++)
        {
          i = TimeBinsStar.ActiveParticleList[idx];

          if(StarNgbs[i] < 1.0)
            {
              npleft++;
              
              SP[i].Hsml *= 2;
            }
          else
            /* Mark as inactive */
            SP[i].DensityFlag = -1; 
        }
     
      sumup_large_ints(1, &npleft, &ntot);

      t1 = second();

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("STAR_DENSITY: ngb iteration %3d: need to repeat for %12lld particles. (took %g sec)\n", iter, ntot,
                       timediff(t0, t1));

          if(iter > MAXITER)
            terminate("failed to converge in neighbour iteration in star_density()\n");
        }
    }
  while(ntot > 0);

  /* 
  Pass 2 — Add star feedback properties to the host cell.
  */

 /* Re-activate all stars */
  for(i = 0; i < NumStars; i++)
    SP[i].DensityFlag = 1;

  pass++;

  generic_comm_pattern(TimeBinsStar.NActiveParticles, kernel_local, kernel_imported);

  myfree(StarHostCell);
  myfree(StarNgbs);

  /* Mark as active */
  for(i = 0; i < NumStars; i++)
     SP[i].DensityFlag = 1;
  
  /* Collect timing information */
  CPU_Step[CPU_INIT] += measure_time();
}

/*! \brief Inner function of the SPH density calculation
 *
 *  This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 *
 *  \param[in] target Index of particle in local data/import buffer.
 *  \param[in] mode Mode in which function is called (local or impored data).
 *  \param[in] threadid ID of local thread.
 *
 *  \return 0
 */
static int star_density_evaluate(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  int ngbs, hosthydrobin = 0; 
  double h, h2, dx, dy, dz, r, r2, wk; 
  MyDouble *pos;

  MyDouble distance = MAX_REAL_NUMBER;
  int index = -1, task = -1;

  data_in local, *target_data;
  data_out out = {0};

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
  h = target_data->Hsml;
  h2 = h * h;

  if(pass == 2)
    {
      index = target_data->HostCell.Index;
      task = target_data->HostCell.Task;
    }

  ngbs = 0;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      i = Thread[threadid].Ngblist[n];

/* compute star->cell position vectors: posSP-posSphP */
      dx = pos[0] - P[i].Pos[0];
      dy = pos[1] - P[i].Pos[1];
      dz = pos[2] - P[i].Pos[2];

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
          /* Pass 1 */
          if(pass == 1)
            {
              ngbs++;

              r = sqrt(r2);
              
              if(r < distance)
                {
                  distance = r;
                  index = i;
                  task = ThisTask;
                }
            }
 
          /* Pass 2 */
          else 
            {          
              if(i == index && ThisTask == task)
                {                
                  hosthydrobin = P[i].TimeBinHydro;

#if defined(TREE_BASED_TIMESTEPS) && defined(SUPERNOVAE)
                  MyDouble time_sn_yr = target_data->TimeSN_yr;
                  MyDouble physical_age_yr = target_data->PhysicalAge_yr;
          
                  if(time_sn_yr < MAX_REAL_NUMBER)
                    {
                      double E_inject_code = 1e51 / 
                      (All.cf_UnitMass_in_g * All.cf_UnitVelocity_in_cm_per_s * All.cf_UnitVelocity_in_cm_per_s);

                      double unew = SphP[i].Utherm + E_inject_code / P[i].Mass;

                      double t_frac = physical_age_yr / time_sn_yr;
                      t_frac = fmin(fmax(t_frac, 0.0), 1.0);

                      double Csn = SphP[i].Csnd + (sqrt(GAMMA * GAMMA_MINUS1 * unew) - SphP[i].Csnd) * t_frac;
          
                      if(Csn > SphP[i].Csn)
                        SphP[i].Csn = Csn;
                    }
#endif

#if defined(WINDS) || defined(SUPERNOVAE)
                  SphP[i].WindsAndSN[SphP[i].Host++] = target_data->WindsAndSN;
#endif
                }
            }
        }
    }

  if(pass == 1)
    {
      out.Ngbs = ngbs;
      out.HostCell.Distance = distance;
      out.HostCell.Index = index;
      out.HostCell.Task = task;
    }

  if(pass == 2)
    out.HostHydroBin = hosthydrobin;

  /* now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

/* \brief Determines if a SP is active in current timestep.
 *
 *  \param[in] n Index of SP in Particle array
 *
 *  \return 1: SP active; 0: SP not active.
 */
int star_density_isactive(int n)
{
  if(SP[n].DensityFlag < 0)
    return 0;

  return 1;
}