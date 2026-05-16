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
  double E_SN = SphP[i].WindsAndSN.SN_EnergyInject;
  double m_ej = SphP[i].WindsAndSN.SN_MassLoss;
  
  double E51 = E_SN * All.cf_UnitEnergy_in_cgs / 1.0e51;
  
  double n_H  = 0.76 * NgbsDensity * All.cf_UnitDensity_in_cgs / PROTONMASS; 

#ifdef METALS
  double Zsol = fmax(NgbsMetallicity / 0.0127, 0.01);
#else
  double Zsol = 1.0;
#endif

  /* Terminal momentum: Kim & Ostriker (2015) */
  double p_term = 3.0e5 /* Msun km/s */
  * pow(E51, 16. / 17.)
  * pow(n_H, -2. / 17.)
  * pow(Zsol, -0.14);
  
  p_term /= (All.cf_UnitMass_in_Msun * All.cf_UnitVelocity_in_cm_per_s / 1.0e5);
  
  /* Boost momentum */
  double fkin = 0.28;
  double a = m_ej * (sum_of_ngb_masses); 
  double b = (m_ej / (2 * fkin * E));
  

  double fboost = fmin((sqrt(b*b + a) - b) / a, p_term / (sqrt(2 * fkin * E_SN * m_ej)));

  double p_SN = fboost * sqrt(2 * fkin * E_SN * m_ej);


  /* energy-conserving limit */
  double p_egy = sqrt(2.0 * m_ej * E_SN);

  double p0 = fmin(p_egy, p_term);

  Usnr = Esnr * (1- fkin * fboost * (fboost * a + 2 * b));

  /* residual thermal energy: what's left after momentum injection */
  *Eth_out = Usnr;

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
      for(int k = 0; k < 3; k++)
        SphP[i].StarMomentumFeed[k] += SphP[i].WindsAndSN.MassLoss * F_HOST * SphP[i].WindsAndSN.StarVelocity[k];

      double sq_vstar = SphP[i].WindsAndSN.StarVelocity[0]*SphP[i].WindsAndSN.StarVelocity[0] 
      + SphP[i].WindsAndSN.StarVelocity[1]*SphP[i].WindsAndSN.StarVelocity[1] 
      + SphP[i].WindsAndSN.StarVelocity[2]*SphP[i].WindsAndSN.StarVelocity[2];

      double sq_vwind = SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss 
      * SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss;

      SphP[i].StarEnergyFeed += 0.5 * SphP[i].WindsAndSN.MassLoss * F_HOST * (sq_vstar + sq_vwind);
#endif 

#ifdef SUPERNOVAE
      SphP[i].StarMassFeed += SphP[i].WindsAndSN.SN_MassLoss * F_HOST;
#ifdef METALS
      SphP[i].StarMetalsFeed += SphP[i].WindsAndSN.SN_MetalsLoss * F_HOST;
#endif
      for(int k = 0; k < 3; k++)
        SphP[i].StarMomentumFeed[k] += SphP[i].WindsAndSN.SN_MassLoss * F_HOST * SphP[i].WindsAndSN.StarVelocity[k];

      double Eth_SN;
      double p0_SN = compute_p0_SN(i, &Eth_SN);

      SphP[i].StarEnergyFeed += Eth_SN * F_HOST;
#endif 

      /* We will loop over the cell faces 4 times */
      int n_faces = 0;
      int dc_list[MAX_FACES];     

      /* Compute weights */
      double nplus[3], nminus[3], fplus[3], fminus[3];
      
      /* Accumulate helper */
      double Splus[3], Sminus[3];
      for(int k = 0; k < 3; k++)
        Splus[k] = Sminus[k] = 0;

      /* Ngbs properties (Start with host) */
      int Ngbs = 1;
      double NgbsMass, NgbsDensity, NgbsMetallicity;     
      NgbsMass = P[i].Mass * F_HOST;
      NgbsDensity = SphP[i].Density * F_HOST;
      NgbsMetallicity = SphP[i].GasMetallicity * F_HOST;

      /* First pass */  
      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("star_feedback: strange connectivity q=%d MaxNvc=%d\n", q, MaxNvc);

          int dp = DC[q].dp_index;
          int vf = DC[q].vf_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            terminate("Particle < 0");

          if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
            particle -= NumGas;

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

              double omega = 0.5 * (1 - 1 / sqrt(1 + 4 * Mesh.VF[vf].area / (M_PI * nn*nn)));
                    
              for(int k = 0; k < 3; k++)
                {
                  nplus[k] = n[k] >= 0 ? n[k] : 0;
                  nminus[k] = n[k] < 0 ? n[k] : 0; 
                        
                  Splus[k] += omega * fabs(nplus[k]);
                  Sminus[k] += omega * fabs(nminus[k]);
                }
            }
              
          if(Mesh.DP[dp].task == ThisTask)
            {
              Ngbs++;
              NgbsMass += SphP[particle].Density * SphP[particle].Volume * (1-F_HOST);
              NgbsDensity += SphP[particle].Density * (1-F_HOST);
              NgbsMetallicity += SphP[particle].GasMetallicity * (1-F_HOST);
            }
          else
            {
              Ngbs++;
              NgbsMass += PrimExch[particle].Density * PrimExch[particle].Volume * (1-F_HOST);
              NgbsDensity += PrimExch[particle].Density * (1-F_HOST);
              NgbsMetallicity += PrimExch[particle].Scalars[METALS_INDEX] * (1-F_HOST);
            }     
          
          if(n_faces >= MAX_FACES)
            terminate("star_feedback: MAX_FACES exceeded for cell %d\n", i);

          dc_list[n_faces++] = q;

          if(q == SphP[i].last_connection) 
            break;
          
          q = DC[q].next;
        }

      NgbsMass /= Ngbs;
      NgbsDensity /= Ngbs;
      NgbsMetallicity /= Ngbs;     

      double weights[MAX_FACES][3];
      double wtot = 0.0; 
      double w[3];
      
      /* Second pass */
      for(int f = 0; f < n_faces; f++)
        {
          q = dc_list[f];
          
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
                
          for(int k = 0; k < 3; k++)
            weights[f][k] = w[k];
                                    
            wtot += sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
        }

      if(wtot <= 0.0)
        terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);
      
      /* Helpers for SN momentum */
      double a, b, num, den = 0;

      /* Third pass */ 
      for(int f = 0; f < n_faces; f++)
        {
          q = dc_list[f];

          int dp = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            terminate("Particle < 0");

          if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
            particle -= NumGas;

          double wbar[3]; 
          
          for(int k = 0; k < 3; k++)
            wbar[k] = weights[f][k] / wtot;

          sq_wbar = (wbar[0]*wbar[0] + wbar[1]*wbar[1] + wbar[2]*wbar[2]);
        
          num = sq_wbar;
          
          if(Mesh.DP[dp].task == ThisTask)
            den = SphP[particle].Density * SphP[particle].Volume + sqrt(sq_wbar)*m_ej*(1-F_HOST);
          else
            den = PrimExch[particle].Density * PrimExch[particle].Volume + sqrt(sq_wbar)*m_ej*(1-F_HOST);
          
          a += num / den;

          if(Mesh.DP[dp].task == ThisTask)
            {
              num = SphP[particle].Density * SphP[particle].Volume 
              * ((P[particle].Vel[0] - SphP[i].WindsAndSN.StarVelocity[0]) * wbar[0] 
              + (P[particle].Vel[1] - SphP[i].WindsAndSN.StarVelocity[1]) * wbar[1]
              + (P[particle].Vel[2] - SphP[i].WindsAndSN.StarVelocity[2]) * wbar[2]);
              den = SphP[particle].Density * SphP[particle].Volume + sqrt(sq_wbar)*m_ej*(1-F_HOST);
            }
          else
            {
              num = PrimExch[particle].Density * PrimExch[particle].Volume 
              * ((PrimExch[particle].VelGas[0] - SphP[i].WindsAndSN.StarVelocity[0]) * wbar[0] 
              + (PrimExch[particle].VelGas[1] - SphP[i].WindsAndSN.StarVelocity[1]) * wbar[1]
              + (PrimExch[particle].VelGas[2] - SphP[i].WindsAndSN.StarVelocity[2]) * wbar[2]);
              den = PrimExch[particle].Density * PrimExch[particle].Volume + sqrt(sq_wbar)*m_ej*(1-F_HOST);
            }

          b += num / den;
        }

      a *= m_ej;

      b *= sqrt(m_ej / (2 * fkin * Esnr));
 
      /* Fourth pass */  
      for(int f = 0; f < n_faces; f++)
        {
          q = dc_list[f];

          double wbar[3]; 
          
          for(int k = 0; k < 3; k++)
            wbar[k] = weights[f][k] / wtot;
      
          /* Outward kick direction: host cell center -> face centroid */
          double n[3], nn; 

          n[0] = Mesh.DP[dp].x - P[i].Pos[0];
          n[1] = Mesh.DP[dp].y - P[i].Pos[1];
          n[2] = Mesh.DP[dp].z - P[i].Pos[2];
          nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

          struct FBKick kick = {0};
          kick.CellIndex = DC[q].index;

#ifdef WINDS
          kick.DeltaMass_wind = SphP[i].WindsAndSN.MassLoss * sqrt(wbar*wbar) * (1-F_HOST);
#ifdef METALS
          kick.DeltaMetals_wind = SphP[i].WindsAndSN.MetalsLoss * sqrt(wbar*wbar) * (1-F_HOST);
#endif    
          for(int k = 0; k < 3; k++)
            kick.DeltaP_wind[k] = SphP[i].WindsAndSN.MassLoss * sqrt(wbar*wbar) * (1-F_HOST) 
            * (SphP[i].WindsAndSN.StarVelocity[k] + SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss * wbar[k]);
          
          double sq_vstar = SphP[i].WindsAndSN.StarVelocity[0]*SphP[i].WindsAndSN.StarVelocity[0] 
          + SphP[i].WindsAndSN.StarVelocity[1]*SphP[i].WindsAndSN.StarVelocity[1] 
          + SphP[i].WindsAndSN.StarVelocity[2]*SphP[i].WindsAndSN.StarVelocity[2];

          double sq_vwind = SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss
          * SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss; 

          double cross = 2 * SphP[i].WindsAndSN.StarVelocity[0] * SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss * wbar[0] 
          + SphP[i].WindsAndSN.StarVelocity[1] * SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss * wbar[1] 
          + SphP[i].WindsAndSN.StarVelocity[2] * SphP[i].WindsAndSN.WindMomentum / SphP[i].WindsAndSN.MassLoss * wbar[2];

          kick.DeltaE_wind = 0.5 * SphP[i].WindsAndSN.MassLoss * sqrt(wbar*wbar) * (1-F_HOST) 
          * (sq_vstar + sq_vwind * wbar*wbar + cross);
          
#endif
 
#ifdef SUPERNOVAE
          kick.DeltaMass_SN = SphP[i].WindsAndSN.SN_MassLoss * sqrt(wbar*wbar) * (1-F_HOST);
#ifdef METALS
          kick.DeltaMetals_SN = SphP[i].WindsAndSN.SN_MetalsLoss * sqrt(wbar*wbar) * (1-F_HOST);
#endif 
          for(int k = 0; k < 3; k++)
            kick.DeltaP_SN[k] = SphP[i].WindsAndSN.SNMassLoss * sqrt(wbar*wbar) * (1-F_HOST) * SphP[i].WindsAndSN.StarVelocity[k]
            + psnr * wbar[k];

          kick.DeltaE_wind = usnr * sqrt(wbar*wbar) * (1-F_HOST);
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
      SphP[i].WindsAndSN[SphP[i].Host--] = {0};
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