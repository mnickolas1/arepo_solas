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

struct Feedback_Kick
{
  int CellIndex; /* cell index on the receiving task */

#ifdef WINDS
  MyDouble DeltaMass;
#ifdef METALS
  MyDouble DeltaMetals;
#endif
  MyDouble DeltaP[3];
  MyDouble DeltaE;
#endif

#ifdef SUPERNOVAE
  MyDouble SN_DeltaMass;
#ifdef METALS
  MyDouble SN_DeltaMetals;
#endif
  MyDouble SN_DeltaP[3];
  MyDouble SN_DeltaE;
#endif
};

/* Apply a kick packet to a local cell */
static void apply_kick(int j, const struct Feedback_Kick *Kick)
{
#ifdef WINDS
      SphP[j].StarMassFeed += Kick->DeltaMass;
#ifdef METALS
      SphP[j].StarMetalsFeed += Kick->DeltaMetals;
#endif
      for(int k = 0; k < 3; k++)
        SphP[j].StarMomentumFeed[k] += Kick->DeltaP[k];

      SphP[j].StarEnergyFeed += Kick->DeltaE;
#endif 

#ifdef SUPERNOVAE
      SphP[j].StarMassFeed += Kick->SN_DeltaMass;
#ifdef METALS
      SphP[j].StarMetalsFeed += Kick->SN_DeltaMetals;
#endif
      for(int k = 0; k < 3; k++)
        SphP[j].StarMomentumFeed[k] += Kick->SN_DeltaP[k];

      SphP[j].StarEnergyFeed += Kick->SN_DeltaE;
#endif 
}

/* Compute SN p0 = min(p_egy, p_term) and residual thermal energy
   i is the index of the host gas cell, which already holds the ambient
   density and metallicity from star_density pass 2. */

#ifdef SUPERNOVAE
static void SN_compute(int i, double NgbsDensity, double NgbsMetallicity, double epsilon, double a, double b, double *p, double *Eth)
{
  Mechanical_Feedback *WindsAndSN = &SphP[i].WindsAndSN[SphP[i].Host - 1];
  
  double E_SN = WindsAndSN->SN_EnergyInject;
  double m_ej = WindsAndSN->SN_MassLoss;

  double sq_vhost = P[i].Vel[0]*P[i].Vel[0] + P[i].Vel[1]*P[i].Vel[1] + P[i].Vel[2]*P[i].Vel[2];
  double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0] 
  + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1] 
  + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];
  double cross = 2.0 * (P[i].Vel[0]*WindsAndSN->StarVelocity[0] 
  + P[i].Vel[1]*WindsAndSN->StarVelocity[1] 
  + P[i].Vel[2]*WindsAndSN->StarVelocity[2]);
  
  double E_SNR = E_SN + 0.5 * (P[i].Mass * WindsAndSN->SN_MassLoss * F_HOST) / (P[i].Mass + WindsAndSN->SN_MassLoss * F_HOST)
  *(sq_vhost + sq_vstar - cross) + 0.5 * epsilon;
  
  double E51 = E_SNR * All.cf_UnitEnergy_in_cgs / 1.0e51;
  
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
  a *= m_ej;
  b *= sqrt(m_ej / (2 * fkin * E_SNR));
  
  double fboost = fmin((sqrt(b*b + a) - b) / a, p_term / (sqrt(2 * fkin * E_SNR * m_ej)));

  double p_SNR = fboost * sqrt(2 * fkin * E_SNR * m_ej);

  double U_SNR = E_SNR * (1 - fkin * fboost * (fboost * a + 2 * b));

  *p = p_SNR;
  *Eth = U_SNR;
}
#endif 

/* 
star_feedback() -> main entry point
*/
void star_feedback(void)
{
  TIMER_START(CPU_STARS_RADIATION);

#define MAX_FACES 128

  int max_export = NumGas * 16;
  struct Feedback_Kick *ExportBuf = mymalloc("ExportBuf",  max_export * sizeof(struct Feedback_Kick));
  int *ExportTask = mymalloc("ExportTask", max_export * sizeof(int));
  int n_export = 0;

  /* Loop over gas cells; act on host cells */
  for(int i = 0; i < NumGas; i++)
    {
      /* Not a host cell -> skip */
      if(SphP[i].Host <= 0) continue; 

      /* We will loop over the cell faces 4 times */
      int n_faces = 0;
      int dc_list[MAX_FACES];     

      /* GEOMETRY: passes 1 & 2 are star-independent (Voronoi mesh) */

      /* Compute weights */
      double nplus[3], nminus[3], fplus[3], fminus[3];
      
      /* Accumulate helper */
      double Splus[3], Sminus[3];
      for(int k = 0; k < 3; k++)
        Splus[k] = Sminus[k] = 0;

      /* First pass */  
      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("star_feedback: strange connectivity q=%d MaxNvc=%d\n", q, MaxNvc);

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

              double omega = 0.5 * (1 - 1 / sqrt(1 + 4 * Mesh.VF[vf].area / (M_PI * nn*nn)));
                    
              for(int k = 0; k < 3; k++)
                {
                  nplus[k] = n[k] >= 0 ? n[k] : 0;
                  nminus[k] = n[k] < 0 ? n[k] : 0; 
                        
                  Splus[k] += omega * fabs(nplus[k]);
                  Sminus[k] += omega * fabs(nminus[k]);
                }
            }
          
          if(n_faces >= MAX_FACES)
            terminate("star_feedback: MAX_FACES exceeded for cell %d\n", i);

          dc_list[n_faces++] = q;

          if(q == SphP[i].last_connection) 
            break;
          
          q = DC[q].next;
        }

      double weights[MAX_FACES][3];
      double wtot = 0.0; 

#ifdef SUPERNOVAE
      /* Ngbs properties (Start with host) */
      int Ngbs = 1;
      double NgbsMass, NgbsDensity, NgbsMetallicity;     
      NgbsMass = P[i].Mass * F_HOST;
      NgbsDensity = SphP[i].Density * F_HOST;
      NgbsMetallicity = SphP[i].GasMetallicity * F_HOST;
#endif
  
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

          double w[3] = {0.0, 0.0, 0.0};

          if(nn > 0.0 && Mesh.VF[vf].area > 0.0)
            {
              n[0] /= nn;  n[1] /= nn;  n[2] /= nn;
                                        
              double omega = 0.5 * (1 - 1 / sqrt(1 + 4 * Mesh.VF[vf].area / (M_PI * nn*nn)));
                    
              for(int k = 0; k < 3; k++)
                {
                  nplus[k] = n[k] >= 0 ? n[k] : 0;
                  nminus[k] = n[k] < 0 ? n[k] : 0; 
                  fplus[k] = (Splus[k]  > 0.0) ? sqrt(0.5 * (1 + Sminus[k]*Sminus[k] / (Splus[k]*Splus[k]))) : 0;
                  fminus[k] =(Sminus[k]  > 0.0) ? sqrt(0.5 * (1 + Splus[k]*Splus[k] / (Sminus[k]*Sminus[k]))) : 0;

                  w[k] = omega * (nplus[k] * fplus[k] + nminus[k] * fminus[k]); 
                }                                         
            }
                
          for(int k = 0; k < 3; k++)
            weights[f][k] = w[k];
                                    
          wtot += sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);

#ifdef SUPERNOVAE
          int particle = Mesh.DP[dp].index;

          if(particle < 0)
            terminate("Particle < 0");

          if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
            particle -= NumGas;

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
#endif    
        }

      if(wtot <= 0.0)
        terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);

#ifdef SUPERNOVAE
      NgbsMass /= Ngbs;
      NgbsDensity /= Ngbs;
      NgbsMetallicity /= Ngbs; 
#endif

      for(int h = 0; h < SphP[i].Host; h++)
        {
          Mechanical_Feedback *WindsAndSN = &SphP[i].WindsAndSN[h];

#ifdef SUPERNOVAE      
          /* Helpers for SN momentum */
          double num, den;
          double epsilon = 0.0, a = 0.0, b = 0.0;
          double m_ej = WindsAndSN->SN_MassLoss;
      
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

              double sq_wbar = (wbar[0]*wbar[0] + wbar[1]*wbar[1] + wbar[2]*wbar[2]);

              double mass, vel[3];
              if(Mesh.DP[dp].task == ThisTask)
                {
                  mass = SphP[particle].Density * SphP[particle].Volume;
                  vel[0] = P[particle].Vel[0];
                  vel[1] = P[particle].Vel[1];
                  vel[2] = P[particle].Vel[2];
                }
              else
                {
                  mass = PrimExch[particle].Density * PrimExch[particle].Volume;
                  vel[0] = PrimExch[particle].VelGas[0];
                  vel[1] = PrimExch[particle].VelGas[1];
                  vel[2] = PrimExch[particle].VelGas[2];
                }

              double sq_vcell = vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2];
              double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0]
              + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1]
              + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];
              double cross  = 2.0 * (vel[0]*WindsAndSN->StarVelocity[0]
              + vel[1]*WindsAndSN->StarVelocity[1]
              + vel[2]*WindsAndSN->StarVelocity[2]);
         
              num = 0.5 * mass * m_ej*sqrt(sq_wbar)*(1-F_HOST) * (sq_vcell + sq_vstar - cross);
              den = mass + m_ej*sqrt(sq_wbar)*(1-F_HOST);
          
              epsilon += num / den;
               
              num = sq_wbar;
              den = mass + m_ej*sqrt(sq_wbar)*(1-F_HOST);
          
              a += num / den;

              num = mass * ((vel[0] - WindsAndSN->StarVelocity[0]) * wbar[0] 
              + (vel[1] - WindsAndSN->StarVelocity[1]) * wbar[1]
              + (vel[2] - WindsAndSN->StarVelocity[2]) * wbar[2]);
              den = mass + m_ej*sqrt(sq_wbar)*(1-F_HOST);

              b += num / den;
            }

          double p, Eth;
          SN_compute(i, NgbsDensity, NgbsMetallicity, epsilon, a, b, &p, &Eth);
#endif
      
          /* Host feedback */
#ifdef WINDS
          SphP[i].StarMassFeed += WindsAndSN->MassLoss * F_HOST;
#ifdef METALS
          SphP[i].StarMetalsFeed += WindsAndSN->MetalsLoss * F_HOST;
#endif
          for(int k = 0; k < 3; k++)
            SphP[i].StarMomentumFeed[k] += WindsAndSN->MassLoss * F_HOST * WindsAndSN->StarVelocity[k];

          double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0] 
          + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1] 
          + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];

          double sq_vwind = WindsAndSN->WindMomentum / WindsAndSN->MassLoss 
          * WindsAndSN->WindMomentum / WindsAndSN->MassLoss;

          SphP[i].StarEnergyFeed += 0.5 * WindsAndSN->MassLoss * F_HOST * (sq_vstar + sq_vwind);
#endif 

#ifdef SUPERNOVAE
          SphP[i].StarMassFeed += WindsAndSN->SN_MassLoss * F_HOST;
#ifdef METALS
          SphP[i].StarMetalsFeed += WindsAndSN->SN_MetalsLoss * F_HOST;
#endif
          for(int k = 0; k < 3; k++)
            SphP[i].StarMomentumFeed[k] += WindsAndSN->SN_MassLoss * F_HOST * WindsAndSN->StarVelocity[k];

          SphP[i].StarEnergyFeed += Eth * F_HOST;
#endif 
 
          /* Fourth pass */  
          for(int f = 0; f < n_faces; f++)
            {
              q = dc_list[f];

              double wbar[3]; 
          
              for(int k = 0; k < 3; k++)
                wbar[k] = weights[f][k] / wtot;

              double sq_wbar = (wbar[0]*wbar[0] + wbar[1]*wbar[1] + wbar[2]*wbar[2]);
      
              struct Feedback_Kick Kick = {0};
              Kick.CellIndex = DC[q].index;
           
              /* Mesh ngbs feedback */
#ifdef WINDS
              Kick.DeltaMass = WindsAndSN->MassLoss * sqrt(sq_wbar) * (1-F_HOST);
#ifdef METALS
              Kick.DeltaMetals = WindsAndSN->MetalsLoss * sqrt(sq_wbar) * (1-F_HOST);
#endif    
              for(int k = 0; k < 3; k++)
                Kick.DeltaP[k] = WindsAndSN->MassLoss * sqrt(sq_wbar) * (1-F_HOST) 
                * (WindsAndSN->StarVelocity[k] + WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[k]);
          
              double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0] 
              + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1] 
              + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];

              double sq_vwind = WindsAndSN->WindMomentum / WindsAndSN->MassLoss
              * WindsAndSN->WindMomentum / WindsAndSN->MassLoss; 

              double cross = 2.0 * (WindsAndSN->StarVelocity[0] * WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[0] 
              + WindsAndSN->StarVelocity[1] * WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[1] 
              + WindsAndSN->StarVelocity[2] * WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[2]);

              Kick.DeltaE = 0.5 * WindsAndSN->MassLoss * sqrt(sq_wbar) * (1-F_HOST) 
              * (sq_vstar + sq_vwind * sq_wbar + cross);   
#endif
 
#ifdef SUPERNOVAE
              Kick.SN_DeltaMass = WindsAndSN->SN_MassLoss * sqrt(sq_wbar) * (1-F_HOST);
#ifdef METALS
              Kick.SN_DeltaMetals = WindsAndSN->SN_MetalsLoss * sqrt(sq_wbar) * (1-F_HOST);
#endif 
              for(int k = 0; k < 3; k++)
                Kick.SN_DeltaP[k] = WindsAndSN->SN_MassLoss * sqrt(sq_wbar) * (1-F_HOST) * WindsAndSN->StarVelocity[k]
                + p * wbar[k];

              Kick.SN_DeltaE = Eth * sqrt(sq_wbar) * (1-F_HOST);
#endif

              if(DC[q].task == ThisTask)
                apply_kick(DC[q].index, &Kick);
              else
                {
                  if(n_export >= max_export)
                    terminate("star_feedback: Feedback_Kick export buffer overflow\n");
                  ExportBuf[n_export] = Kick;
                  ExportTask[n_export] = DC[q].task;
                  n_export++;
                }
            }

          /* Clear this star's feedback entry */
          memset(WindsAndSN, 0, sizeof(*WindsAndSN));
        }
      
      /* All stars processed: release host slot */        
      SphP[i].Host = 0;
    }

  /* MPI exchange of remote kick packets via MPI_Alltoallv */
  int *SendCount = mymalloc("FBSendCount", NTask * sizeof(int));
  int *RecvCount = mymalloc("FBRecvCount", NTask * sizeof(int));
  int *SendDisp = mymalloc("FBSendDisp",  NTask * sizeof(int));
  int *RecvDisp = mymalloc("FBRecvDisp",  NTask * sizeof(int));
 
  memset(SendCount, 0, NTask * sizeof(int));
  for(int k = 0; k < n_export; k++)
    SendCount[ExportTask[k]]++;
 
  MPI_Alltoall(SendCount, 1, MPI_INT, RecvCount, 1, MPI_INT, MPI_COMM_WORLD);
 
  SendDisp[0] = RecvDisp[0] = 0;
  for(int t = 1; t < NTask; t++)
    {
      SendDisp[t] = SendDisp[t-1] + SendCount[t-1];
      RecvDisp[t] = RecvDisp[t-1] + RecvCount[t-1];
    }
  int n_recv = RecvDisp[NTask-1] + RecvCount[NTask-1];
 
  /* Sort ExportBuf into task-contiguous order for Alltoallv */
  struct Feedback_Kick *SortedExport = mymalloc("FBSortedExport", (n_export > 0 ? n_export : 1) * sizeof(struct Feedback_Kick));
  
  int *tmp_offset = mymalloc("FBTmpOffset", NTask * sizeof(int));
  memcpy(tmp_offset, SendDisp, NTask * sizeof(int));
  for(int k = 0; k < n_export; k++)
    {
      int t = ExportTask[k];
      SortedExport[tmp_offset[t]++] = ExportBuf[k];
    }
  myfree(tmp_offset);
 
  struct Feedback_Kick *RecvBuf = mymalloc("RecvBuf",
  (n_recv > 0 ? n_recv : 1) * sizeof(struct Feedback_Kick));
 
  const int sz = (int)sizeof(struct Feedback_Kick);
 
  int *SendCountB = mymalloc("SendCountB", NTask * sizeof(int));
  int *RecvCountB = mymalloc("RecvCountB", NTask * sizeof(int));
  int *SendDispB  = mymalloc("SendDispB",  NTask * sizeof(int));
  int *RecvDispB  = mymalloc("RecvDispB",  NTask * sizeof(int));
 
  for(int t = 0; t < NTask; t++)
    {
      SendCountB[t] = SendCount[t] * sz;
      RecvCountB[t] = RecvCount[t] * sz;
      SendDispB[t]  = SendDisp[t]  * sz;
      RecvDispB[t]  = RecvDisp[t]  * sz;
    }
 
  MPI_Alltoallv(SortedExport, SendCountB, SendDispB, MPI_BYTE,
  RecvBuf, RecvCountB, RecvDispB, MPI_BYTE,
  MPI_COMM_WORLD);
 
  /* Free byte-count arrays in reverse allocation order (LIFO stack) */
  myfree(RecvDispB);
  myfree(SendDispB);
  myfree(RecvCountB);
  myfree(SendCountB);
 
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
 
  TIMER_STOP(CPU_STARS_FEEDBACK);
}