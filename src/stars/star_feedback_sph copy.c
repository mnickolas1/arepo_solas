#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

int star_mesh_export(int target, int task, int thread_id)
{
  if(Thread[thread_id].Exportflag[task] != target)
    {
      Thread[thread_id].Exportflag[task] = target;
      int nexp = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task  = task;
      Thread[thread_id].PartList[nexp].Index = target;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp = Thread[thread_id].NexportNodes++;
  nexp = -1 - nexp;
  struct data_partlist *partlist = (struct data_partlist *)(((char *)Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task = task;
  nodelist[nexp].Index = target;

  Thread[thread_id].ExportSpace -= sizeof(struct data_partlist) + sizeof(int);
  return 0;
}


void star_apply_feedback(void)
{
#ifdef WINDS
      SphP[i].StarMassFeed += SphP[i].WindsAndSN.MassLoss * F_HOST;
#ifdef METALS
      SphP[i].StarMetalsFeed += SphP[i].WindsAndSN.MetalsLoss * F_HOST;
#endif
      for(int k = 0; k < 3; k++)
        SphP[i].StarMomentumFeed[k] += SphP[i].WindsAndSN.MassLoss * F_HOST * SphP[i].WindsAndSN.StarVelocity[k];
#endif 

#ifdef SUPERNOVAE
      SphP[i].StarMassFeed += SphP[i].WindsAndSN.SN_MassLoss * F_HOST;
#ifdef METALS
      SphP[i].StarMetalsFeed += SphP[i].WindsAndSN.SN_MetalsLoss * F_HOST;
#endif
      double Eth_SN;
      double p0_SN = compute_p0_SN(i, &Eth_SN);

      SphP[i].StarEnergyFeed += Eth_SN * F_HOST;
#endif 

#ifdef WINDS
          kick.DeltaMass_wind = SphP[i].WindsAndSN.MassLoss * sqrt(wbar*wbar) * (1 - F_HOST);
#ifdef METALS
          kick.DeltaMetals_wind = SphP[i].WindsAndSN.MetalsLoss * sqrt(wbar*wbar) * (1 - F_HOST);
#endif    
          for(int k = 0; k < 3; k++)
            kick.DeltaP_wind[0] = SphP[i].WindsAndSN.WindMomentum * wbar;
#endif
 
#ifdef SUPERNOVAE
          kick.DeltaMass_SN = SphP[i].WindsAndSN.SN_MassLoss * sqrt(wbar*wbar) * (1 - F_HOST);
#ifdef METALS
          kick.DeltaMetals_SN = SphP[i].WindsAndSN.SN_MetalsLoss * sqrt(wbar*wbar) * (1 - F_HOST);
#endif 
          kick.DeltaP_SN[0] = wbar * p0_SN * nx;
          kick.DeltaP_SN[1] = wbar * p0_SN * ny;
          kick.DeltaP_SN[2] = wbar * p0_SN * nz;
#endif
}

void star_mesh_loop(void)
{
  /* Loop over gas cells; act on host cells */
  for(int i = 0; i < NumGas; i++)
    {
      /* Not a host cell -> skip */
      if(SphP[i].Host <= 0) continue;   

      /* Compute weights */
      double nplus[3], nminus[3], fplus[3], fminus[3];
      
      /* Accumulate */
      double Splus[3], Sminus[3];
      for(int k = 0; k < 3; k++)
        Splus[k] = Sminus[k] = 0;

      /* First pass */  
      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("star_feedback: strange connectivity q=%d MaxNvc=%d\n", q, MaxNvc);

            if(DC[q].task >= 0 && DC[q].task < NTask)
              {
                int dp = DC[q].dp_index;
                int vf = DC[q].vf_index;
               
                /* Outward face normal (Voronoi face normal) */
                
                double n[3], nn; 

                n[0] = Mesh.DP[dp].x - P[i].Pos[0];
                n[1] = Mesh.DP[dp].y - P[i].Pos[1];
                n[2] = Mesh.DP[dp].z - P[i].Pos[2];
                nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

                if(nn > 0.0 && Mesh.VF[vf].area > 0.0)
                  {
                    n[0] /= nn;  n[1] /= nn;  n[2] /= nn;

                    double omega = 0.5 * (1 - 1 / sqrt(1 + 4 * Mesh.VF[vf].area / (M_PI * nn*nn)));
                    
                    for(int k = 0; k < 3; k++)
                      {
                        nplus[k] = n[k] >= 0 ? n[k] : 0;
                        nminus[k] = n[k] < 0 ? n[k] : 0; 
                        
                        Splus[k] += omega * fabs(nplus[k]);
                        Sminus[k] += omega * fabs(nminus[k]);
                      }
                  }
              }

          if(q == SphP[i].last_connection) 
            break;
          
          q = DC[q].next;
        }

      int dc_list[MAX_FACES];
      double weights[MAX_FACES][3];
      int n_faces = 0;
      double wtot = 0.0; 
      
      double w[3];
      
      /* Second pass */
      q = SphP[i].first_connection;
      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("star_feedback: strange connectivity q=%d MaxNvc=%d\n", q, MaxNvc);

            if(DC[q].task >= 0 && DC[q].task < NTask)
              {
                int dp = DC[q].dp_index;
                int vf = DC[q].vf_index;

                double n[3], nn; 

                n[0] = Mesh.DP[dp].x - P[i].Pos[0];
                n[1] = Mesh.DP[dp].y - P[i].Pos[1];
                n[2] = Mesh.DP[dp].z - P[i].Pos[2];
                nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

                if(nn > 0.0 && Mesh.VF[vf].area > 0.0)
                  {
                    n[0] /= nn;  n[1] /= nn;  n[2] /= nn;
                                        
                    double omega = 0.5 * (1 - 1 / sqrt(1 + 4 * Mesh.VF[vf].area / (M_PI * nn*nn)));
                    
                    for(int k = 0; k < 3; k++)
                      {
                        nplus[k] = n[k] >= 0 ? n[k] : 0;
                        nminus[k] = n[k] < 0 ? n[k] : 0; 
                        fplus[k] = sqrt(0.5 * (1 + Sminus[k]*Sminus[k] / (Splus[k]*Splus[k])));
                        fminus[k] = sqrt(0.5 * (1 + Splus[k]*Splus[k] / (Sminus[k]*Sminus[k])));

                        w[k] = omega * (nplus[k] * fplus[k] + nminus[k] * fminus[k]); 
                      }                                         
                  }

                if(n_faces >= MAX_FACES)
                  terminate("star_feedback: MAX_FACES exceeded for cell %d\n", i);

                dc_list[n_faces] = q;
                
                for(int k = 0; k < 3; k++)
                  weights[n_faces][k] = w[k];
                                      
                wtot += sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

                n_faces++;
              }

          if(q == SphP[i].last_connection) 
            break;
          
          q = DC[q].next;
        }

      if(wtot <= 0.0)
        terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);
      
      /* Third pass */ 
      for(int f = 0; f < n_faces; f++)
        {
          q = dc_list[f];

          double wbar[3]; 
          
          for(int k = 0; k < 3; k++)
            wbar[k] = weights[f][k] / wtot;
          
          if(wbar <= 0.0) 
            terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);

          /* Outward kick direction: host cell center -> face centroid */
          double n[3], nn; 

          n[0] = Mesh.DP[dp].x - P[i].Pos[0];
          n[1] = Mesh.DP[dp].y - P[i].Pos[1];
          n[2] = Mesh.DP[dp].z - P[i].Pos[2];
          nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        }
      /* CLEAR - reset host-cell feedback flags so they don't fire again */
      SphP[i].WindsAndSN[SphP[i].Host--] = {0};
    }
}

static int star_feedback_evaluate(int target, int mode, int threadid);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];

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
  int i, j;

  int threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      i = NextParticle++;

      if(i >= NumGas)
        break;
      
      if(SphP[i].Host > 0)    
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

  generic_comm_pattern(NumGas, kernel_local, kernel_imported);

  star_mesh_loop();
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

  
  if(mode == 0)
    {
      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("star_feedback: strange connectivity q=%d MaxNvc=%d\n", q, MaxNvc);

          int dp = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;

          if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
            particle -= NumGas;

          if(DC[q].task == ThisTask)
            {
              ngbs++;
              ngbsdensity += SphP[particle].Density;
              ngbsmetallicity += SphP[particle].Metallicity;
            }
          
          if(DC[q].task != ThisTask)
            star_mesh_export(i, DC[q].task, threadid);
          
        }
    }
  
  if(mode == 1)
    

  return 0;
}