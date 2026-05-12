#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* Fraction of feedback scalars that go to the host cell */
#define F_HOST (1.0 / 20.0) 

/* 
   star_feedback.c

   Loops over gas cells. Any cell with SphP[i].Host > 0 was marked by
   star_density() as hosting a star, and already carries the feedback
   quantities (MassLoss, WindMomentum, SN_EnergyInject, ...).
*/

/* Kick packet sent to remote face-neighbor cells
   Carries momentum contributions from each active feedback channel. */

struct FBKick
{
  int CellIndex; /* cell index on the receiving task */

#ifdef WINDS
  MyDouble DeltaP_wind[3];
#endif

#ifdef SUPERNOVAE
  MyDouble DeltaP_SN[3];
#endif
};

/* Apply a kick packet to a local cell */
static void apply_kick(int j, const struct FBKick *k)
{
#ifdef WINDS
  P[j].Vel[0] += k->DeltaP_wind[0] / P[j].Mass;
  P[j].Vel[1] += k->DeltaP_wind[1] / P[j].Mass;
  P[j].Vel[2] += k->DeltaP_wind[2] / P[j].Mass;
#endif

#ifdef SUPERNOVAE
  P[j].Vel[0] += k->DeltaP_SN[0] / P[j].Mass;
  P[j].Vel[1] += k->DeltaP_SN[1] / P[j].Mass;
  P[j].Vel[2] += k->DeltaP_SN[2] / P[j].Mass;
#endif
}

/* Compute SN p0 = min(p_egy, p_term) and residual thermal energy
   i is the index of the host gas cell, which already holds the ambient
   density and metallicity from star_density pass 2. */

#ifdef SUPERNOVAE
static double compute_p0_SN(int i, double *Eth_out)
{
  double E_SN = SphP[i].SN_EnergyInject;
  double m_ej = SphP[i].SN_MassLoss;

  double n_H  = SphP[i].Density * All.cf_UnitDensity_in_cgs / PROTONMASS;
#ifdef METALS
  double Zsol = fmax(SphP[i].Metallicity / SOLAR_METALLICITY, 0.01);
#else
  double Zsol = 1.0;
#endif

  /* terminal momentum: Kim & Ostriker (2015) */
  double E51 = E_SN * All.cf_UnitEnergy_in_cgs / 1.0e51;
  double p_term = 2.8e5  /* Msun km/s */
                  * pow(E51,                0.93)
                  * pow(fmax(n_H, 1.0e-4), -0.17)
                  * pow(Zsol,              -0.14);
  p_term /= (All.cf_UnitMass_in_Msun * All.cf_UnitVelocity_in_cm_per_s / 1.0e5);

  /* energy-conserving limit */
  double p_egy = sqrt(2.0 * m_ej * E_SN);

  double p0 = fmin(p_egy, p_term);

  /* residual thermal energy: what's left after momentum injection */
  *Eth_out = fmax(0.0, E_SN - (p0 * p0) / (2.0 * P[i].Mass));

  return p0;
}
#endif 

/* 
star_feedback() -> main entry point
*/
void star_feedback(void)
{
  CPU_Step[CPU_MISC] += measure_time();

#define MAX_FACES 128

  int max_export = NumGas * 16;
  struct FBKick *ExportBuf = mymalloc("FBExportBuf",  max_export * sizeof(struct FBKick));
  int *ExportTask = mymalloc("FBExportTask", max_export * sizeof(int));
  int n_export = 0;

  /* Loop over gas cells; act on host cells */
  for(int i = 0; i < NumGas; i++)
    {
      /* Not a host cell -> skip */
      if(SphP[i].Host <= 0) continue;   

#ifdef WINDS
      SphP[i].StarMassFeed += SphP[i].WindsAndSN.MassLoss * F_HOST;
#ifdef METALS
      SphP[i].StarMetalsFeed += SphP[i].WindsAndSN.MetalsLoss * F_HOST;
#endif
      SphP[i].StarMomentumFeed += SphP[i].WindsAndSN.MassLoss * F_HOST * Vstar;
#endif 

#ifdef SUPERNOVAE
      SphP[i].StarMassFeed += SphP[i].WindsAndSN.SN_MassLoss * F_HOST;
#ifdef METALS
      SphP[i].StarMetalsFeed += SphP[i].WindsAndSN.SN_MetalsLoss * F_HOST;
#endif
      double Eth_SN;
      double p0_SN = compute_p0_SN(i, &Eth_SN);

      SphP[i].Energy += Eth_SN * F_HOST;
#endif 

      /* Compute weights */
      double nplus[3], nminus[3], Splus[3], Sminus[3], fplus[3], fminus[3];
      
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
               
                /* Outward face normal: vector from cell generator to neighbor
                generator, normalized (this is the Voronoi face normal) */
                
                double n[3], nn; 

                n[0] = Mesh.DP[dp].x - P[i].Pos[0];
                n[1] = Mesh.DP[dp].y - P[i].Pos[1];
                n[2] = Mesh.DP[dp].z - P[i].Pos[2];
                nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

                if(nn > 0.0 && Mesh.VF[vf].area > 0.0)
                  {
                    n[0] /= nn;  n[1] /= nn;  n[2] /= nn;

                    double omega = 0.5 * (1 - 1 / sqrt(1 + 4 * Mesh.VF[vf].area / (M_PI * nn*nn)))
                    
                    for(int k = 0; k < 3; k++)
                      {
                        nplus[k] = n[k] >= 0 ? n[k] : 0;
                        nminus[k] = n[k] < 0 ? n[k] : 0; 
                        
                        Splus[k] += omega * nplus[k];
                        Sminus[k] += omega * nminus[k];
                      }
                  }
              }

          if(q == SphP[i].last_connection) 
            break;
          
          q = DC[q].next;
        }

      int dc_list[MAX_FACES];
      double weights[MAX_FACES];
      int n_faces = 0;
      double wtot = 0.0;    
      
      /* Second pass */
      int q = SphP[i].first_connection;
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
                
                double w = 0;

                if(nn > 0.0 && Mesh.VF[vf].area > 0.0)
                  {
                    n[0] /= nn;  n[1] /= nn;  n[2] /= nn;
                                        
                    double omega = 0.5 * (1 - 1 / sqrt(1 + 4 * Mesh.VF[vf].area / (M_PI * nn*nn)));
                    
                    for(int k = 0; k < 3; k++)
                      {
                        fplus[k] = sqrt(0.5 * (1 + Sminus*Sminus / (Splus*Splus)));
                        fminus[k] = sqrt(0.5 * (1 + Splus*Splus / (Sminus*Sminus)));

                        w += omega * (nplus[k] * fplus[k] + nminus[k] * fminus[k]); 
                      }                                         
                  }

                if(n_faces >= MAX_FACES)
                  terminate("star_feedback: MAX_FACES exceeded for cell %d\n", i);

                dc_list[n_faces] = q;
                weights[n_faces] = w;
                wtot += w;
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

          double wbar = weights[f] / wtot;
          
          if(wbar <= 0.0) 
            terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);

          /* Outward kick direction: host cell center -> face centroid */
          double dx = Mesh.VF[vf].cx - P[i].Pos[0];
          double dy = Mesh.VF[vf].cy - P[i].Pos[1];
          double dz = Mesh.VF[vf].cz - P[i].Pos[2];

          double r = sqrt(dx*dx + dy*dy + dz*dz);

          struct FBKick kick = {0};
          kick.CellIndex = DC[q].index;

#ifdef WINDS
          kick.DeltaP_wind[0] = wbar * SphP[i].WindMomentum * nx;
          kick.DeltaP_wind[1] = wbar * SphP[i].WindMomentum * ny;
          kick.DeltaP_wind[2] = wbar * SphP[i].WindMomentum * nz;
          kick.DeltaMass_wind = wbar * SphP[i].MassLoss * (1 - F_HOST);
#ifdef METALS
          kick.DeltaMetals_wind = wbar * SphP[i].MetalsLoss * (1 - F_HOST);
#endif
#endif
 
#ifdef SUPERNOVAE
          kick.DeltaP_SN[0] = wbar * p0_SN * nx;
          kick.DeltaP_SN[1] = wbar * p0_SN * ny;
          kick.DeltaP_SN[2] = wbar * p0_SN * nz;
          kick.DeltaMass_SN = wbar * SphP[i].SN_MassLoss * (1 - F_HOST);
#ifdef METALS
          kick.DeltaMetals_SN = wbar * SphP[i].SN_MetalsLoss * (1 - F_HOST);
#endif 
#endif

          if(DC[q].task == ThisTask)
            apply_kick(DC[q].index, &kick);
          else
            {
              if(n_export >= max_export)
                terminate("star_feedback: FBKick export buffer overflow\n");
              ExportBuf[n_export] = kick;
              ExportTask[n_export] = DC[q].task;
              n_export++;
            }
        }

      /* CLEAR - reset host-cell feedback flags so they don't fire again */

      SphP[i].Host = 0;
      SphP[i].WindsAndSN = {0};
    }

  /* Exchange remote kick packets via MPI_Alltoallv */
  int *SendCount = mymalloc("FBSendCount", NTask * sizeof(int));
  int *RecvCount = mymalloc("FBRecvCount", NTask * sizeof(int));
  int *SendDisp = mymalloc("FBSendDisp",  NTask * sizeof(int));
  int *RecvDisp = mymalloc("FBRecvDisp",  NTask * sizeof(int));

  memset(SendCount, 0, NTask * sizeof(int));
  for(int k = 0; k < n_export; k++)
    SendCount[ExportTask[k]]++;

  MPI_Alltoall(SendCount, 1, MPI_INT,
               RecvCount, 1, MPI_INT, MPI_COMM_WORLD);

  SendDisp[0] = RecvDisp[0] = 0;
  for(int t = 1; t < NTask; t++)
    {
      SendDisp[t] = SendDisp[t-1] + SendCount[t-1];
      RecvDisp[t] = RecvDisp[t-1] + RecvCount[t-1];
    }
  int n_recv = RecvDisp[NTask-1] + RecvCount[NTask-1];

  /* Sort ExportBuf into task-contiguous order */
  struct FBKick *SortedExport = mymalloc("FBSortedExport", fmax(n_export, 1) * sizeof(struct FBKick));
  {
    int *tmp_offset = mymalloc("FBTmpOffset", NTask * sizeof(int));
    memcpy(tmp_offset, SendDisp, NTask * sizeof(int));
    for(int k = 0; k < n_export; k++)
      {
        int t = ExportTask[k];
        SortedExport[tmp_offset[t]++] = ExportBuf[k];
      }
    myfree(tmp_offset);
  }

  struct FBKick *RecvBuf = mymalloc("FBRecvBuf", fmax(n_recv, 1) * sizeof(struct FBKick));

  MPI_Alltoallv(SortedExport, SendCount, SendDisp, MPI_BYTE,
                RecvBuf,      RecvCount, RecvDisp, MPI_BYTE,
                MPI_COMM_WORLD);

  /* Apply received kicks to local cells */
  for(int k = 0; k < n_recv; k++)
    apply_kick(RecvBuf[k].CellIndex, &RecvBuf[k]);

  /* Cleanup in reverse allocation order */
  myfree(RecvBuf);
  myfree(SortedExport);
  myfree(RecvDisp);
  myfree(SendDisp);
  myfree(RecvCount);
  myfree(SendCount);
  myfree(ExportTask);
  myfree(ExportBuf);

  CPU_Step[CPU_FEEDBACK] += measure_time();
}