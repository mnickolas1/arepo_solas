#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"


/* Fraction of feedback scalars that go to the host cell */
#define F_HOST (1.0 / 20.0) 

/* Loops over gas cells. Any cell with SphP[i].Host > 0 was marked by
   star_density() as hosting a star, and already carries the feedback
   quantities (MassLoss, WindMomentum, SN_EnergyInject, ...) */

/* Kick packet sent to remote face-neighbor cells
   -carries momentum contributions from each active feedback channel */

struct Feedback_Kick
{ 
  /* Cell index on the receiving task */
  int CellIndex; 

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
  All.StarFeedbackLocal[0] += Kick->DeltaMass;
#ifdef METALS
  SphP[j].StarMetalsFeed += Kick->DeltaMetals;
  All.StarFeedbackLocal[1] += Kick->DeltaMetals;
#endif
  for(int k = 0; k < 3; k++)
    SphP[j].StarMomentumFeed[k] += Kick->DeltaP[k];

  SphP[j].StarEnergyFeed += Kick->DeltaE;
  All.StarFeedbackLocal[2] += Kick->DeltaE;
#endif 

#ifdef SUPERNOVAE
  SphP[j].StarMassFeed += Kick->SN_DeltaMass;
  All.StarFeedbackLocal[0] += Kick->SN_DeltaMass;
#ifdef METALS
  SphP[j].StarMetalsFeed += Kick->SN_DeltaMetals;
  All.StarFeedbackLocal[1] += Kick->SN_DeltaMetals;
#endif
  for(int k = 0; k < 3; k++)
    SphP[j].StarMomentumFeed[k] += Kick->SN_DeltaP[k];

  SphP[j].StarEnergyFeed += Kick->SN_DeltaE;
  All.StarFeedbackLocal[2] += Kick->SN_DeltaE;
#endif 
}

#ifdef SUPERNOVAE
/* Compute SN p0 = min(p_egy, p_term) and residual thermal energy
   i is the index of the host gas cell, which already holds the ambient
   density and metallicity from star_density pass 2 */
static void SN_compute(int ev, int h, double e, double a, double b, double NgbsDensity, double NgbsMetallicity, double *p, double *Eth)
{
  int i = MechanicalFeedbackEvents.Data[ev].HostIndex;

  Mechanical_Feedback_Data *data = &MechanicalFeedbackEvents.Data[ev + h];
  Mechanical_Feedback *WindsAndSN = &data->WindsAndSN;
  
  double E_SN = WindsAndSN->SN_EnergyInject;
  double m_ej = WindsAndSN->SN_MassLoss;

  /* Pure thermal injection: no ejecta mass, skip momentum prescription */
  if(m_ej <= 0.0)
    {
      *p = 0.0;
      *Eth = E_SN; 
      
      return;
    }

  double mhost, vhost[3];

  mhost = P[i].Mass + SphP[i].StarMassFeed;
  
  for(int k = 0; k < 3; k++)
    vhost[k] = (SphP[i].Momentum[k] + SphP[i].StarMomentumFeed[k]) / mhost;

  double sq_vhost = vhost[0]*vhost[0] + vhost[1]*vhost[1] + vhost[2]*vhost[2];

  double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0] 
  + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1] 
  + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];

  double cross = 2.0 * (vhost[0]*WindsAndSN->StarVelocity[0] 
  + vhost[1]*WindsAndSN->StarVelocity[1] 
  + vhost[2]*WindsAndSN->StarVelocity[2]);
  
  double E_SNR = E_SN + 0.5 * (mhost * WindsAndSN->SN_MassLoss * F_HOST) / (mhost + WindsAndSN->SN_MassLoss * F_HOST)
  * (sq_vhost + sq_vstar - cross) + e;
  
  double E51 = E_SNR * All.cf_UnitEnergy_in_cgs / 1.0e51;
  
  double n_H  = 0.76 * NgbsDensity * All.cf_UnitDensity_in_cgs / PROTONMASS; 

  double Zsol = fmax(NgbsMetallicity / 0.0127, 0.01);

  /* Terminal momentum: Kim & Ostriker (2015) */
  double p_term = 3.0e5 /* Msun km/s */
  * pow(E51, 16.0 / 17.0)
  * pow(n_H, -2.0 / 17.0)
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

/* star_feedback() -> main entry point */
void star_feedback(void)
{
  TIMER_START(CPU_STARS_FEEDBACK);

  #define MAX_FACES 128

  int ev, h, i, k, q, f;

  int n_export = 0;
  int max_export = 100 * NTask;
  int *ExportTask = mymalloc("ExportTask", max_export * sizeof(int));
  struct Feedback_Kick *ExportBuf = mymalloc("ExportBuf",  max_export * sizeof(struct Feedback_Kick));

  /* Act on host cells */
  for(ev = 0; ev < MechanicalFeedbackEvents.NumEvents;)
    {
      i = MechanicalFeedbackEvents.Data[ev].HostIndex;

      int flag_winds_any = 0, flag_sn_any = 0;
      for(h = 0; h < SphP[i].Host; h++)
        {
          Mechanical_Feedback_Data *data = &MechanicalFeedbackEvents.Data[ev + h];
          Mechanical_Feedback *WindsAndSN = &data->WindsAndSN;

#ifdef WINDS
          if(WindsAndSN->MassLoss)
            flag_winds_any++;
#endif

#ifdef SUPERNOVAE
          if(WindsAndSN->SN_MassLoss || WindsAndSN->SN_EnergyInject)
            flag_sn_any++;
#endif
        }
      
      /* This cell only contains dead stars */  
      if(!flag_winds_any && !flag_sn_any)
        {
          /* Go to next host */
          ev += SphP[i].Host;
          /* All stars processed: release host slot */        
          SphP[i].Host = 0;

          continue;
        }

      for(h = 0; h < SphP[i].Host; h++)
        {                    
          int flag_winds = 0, flag_sn = 0;

          Mechanical_Feedback_Data *data = &MechanicalFeedbackEvents.Data[ev + h];
          Mechanical_Feedback *WindsAndSN = &data->WindsAndSN;

#ifdef WINDS
          if(WindsAndSN->MassLoss)
            flag_winds++;
#endif

#ifdef SUPERNOVAE
          if(WindsAndSN->SN_MassLoss || WindsAndSN->SN_EnergyInject)
            flag_sn++;
#endif
          
          /* Dead star */
          if(!flag_winds && !flag_sn)
            continue;

          /* We will loop over the cell faces 4 times */
          int n_faces = 0;
          int dc_list[MAX_FACES];     

          /* GEOMETRY: passes 1 & 2 (Voronoi mesh) */

          /* Compute weights */
          double nplus[3], nminus[3], fplus[3], fminus[3];
      
          /* Accumulate helper */
          double Splus[3], Sminus[3];
          for(k = 0; k < 3; k++)
            Splus[k] = Sminus[k] = 0;

          /* First pass */            
          q = SphP[i].first_connection;
          while(q >= 0)
            {
              int dp = DC[q].dp_index;
              int vf = DC[q].vf_index;
              int particle = Mesh.DP[dp].index;
          
              /* Cell has been removed */
              if(particle < 0)
                {
                  /* Move to next ngb */
                  q = DC[q].next;
                  continue;
                }

              /* Face normal - from cell generator to cell generator */
              double n[3], nn; 

              n[0] = Mesh.DP[dp].x - P[i].Pos[0];
              n[1] = Mesh.DP[dp].y - P[i].Pos[1];
              n[2] = Mesh.DP[dp].z - P[i].Pos[2];
          
              nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

              /* Star-to-face direction */
              double d[3], dd; 
                  
              d[0] = Mesh.VF[vf].cx - WindsAndSN->StarPosition[0];
              d[1] = Mesh.VF[vf].cy - WindsAndSN->StarPosition[1];
              d[2] = Mesh.VF[vf].cz - WindsAndSN->StarPosition[2];
              
              dd = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

              if(nn > 0.0 && dd > 0.0 && Mesh.VF[vf].area > 0.0)
                {
                  n[0] /= nn;  n[1] /= nn;  n[2] /= nn;
                  
                  d[0] /= dd;  d[1] /= dd;  d[2] /= dd;

                  double costheta = n[0]*d[0] + n[1]*d[1] + n[2]*d[2];
                  if(costheta < 0.0) costheta = 0.0;
            
                  double omega = 0.5 * (1 - 1 / sqrt(1 + Mesh.VF[vf].area * costheta / (M_PI * dd*dd)));

                  /* Star-to-cell direction */
                  double r[3], rr; 
                  
                  r[0] = Mesh.DP[dp].x - WindsAndSN->StarPosition[0];
                  r[1] = Mesh.DP[dp].y - WindsAndSN->StarPosition[1];
                  r[2] = Mesh.DP[dp].z - WindsAndSN->StarPosition[2];

                  rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

                  if(rr > 0.0)
                    {
                      r[0] /= rr;  r[1] /= rr;  r[2] /= rr;

                      for(k = 0; k < 3; k++)
                        {
                          nplus[k] = r[k] >= 0 ? r[k] : 0;
                          nminus[k] = r[k] < 0 ? r[k] : 0; 
                        
                          Splus[k] += omega * fabs(nplus[k]);
                          Sminus[k] += omega * fabs(nminus[k]);
                        }
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
  
          /* Second pass */
          for(f = 0; f < n_faces; f++)
            {
              q = dc_list[f];
          
              int dp = DC[q].dp_index;
              int vf = DC[q].vf_index;

              double w[3] = {0.0, 0.0, 0.0};

             /* Face normal - from cell generator to cell generator */
              double n[3], nn; 

              n[0] = Mesh.DP[dp].x - P[i].Pos[0];
              n[1] = Mesh.DP[dp].y - P[i].Pos[1];
              n[2] = Mesh.DP[dp].z - P[i].Pos[2];
          
              nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

              /* Star-to-face direction */
              double d[3], dd; 
                  
              d[0] = Mesh.VF[vf].cx - WindsAndSN->StarPosition[0];
              d[1] = Mesh.VF[vf].cy - WindsAndSN->StarPosition[1];
              d[2] = Mesh.VF[vf].cz - WindsAndSN->StarPosition[2];
              
              dd = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

              if(nn > 0.0 && dd > 0.0 && Mesh.VF[vf].area > 0.0)
                {
                  n[0] /= nn;  n[1] /= nn;  n[2] /= nn;
                  
                  d[0] /= dd;  d[1] /= dd;  d[2] /= dd;

                  double costheta = n[0]*d[0] + n[1]*d[1] + n[2]*d[2];
                  if(costheta < 0.0) costheta = 0.0;
            
                  double omega = 0.5 * (1 - 1 / sqrt(1 + Mesh.VF[vf].area * costheta / (M_PI * dd*dd)));

                  /* Star-to-cell direction */
                  double r[3], rr; 
                  
                  r[0] = Mesh.DP[dp].x - WindsAndSN->StarPosition[0];
                  r[1] = Mesh.DP[dp].y - WindsAndSN->StarPosition[1];
                  r[2] = Mesh.DP[dp].z - WindsAndSN->StarPosition[2];

                  rr = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

                  if(rr > 0.0)
                    {
                      r[0] /= rr;  r[1] /= rr;  r[2] /= rr;
                    
                      for(k = 0; k < 3; k++)
                        {
                          nplus[k] = r[k] >= 0 ? r[k] : 0;
                          nminus[k] = r[k] < 0 ? r[k] : 0; 
                          fplus[k] = (Splus[k]  > 0.0) ? sqrt(0.5 * (1 + Sminus[k]*Sminus[k] / (Splus[k]*Splus[k]))) : 0;
                          fminus[k] =(Sminus[k]  > 0.0) ? sqrt(0.5 * (1 + Splus[k]*Splus[k] / (Sminus[k]*Sminus[k]))) : 0;

                          w[k] = omega * (nplus[k] * fplus[k] + nminus[k] * fminus[k]); 
                        }
                    }                                         
                }
                
              for(k = 0; k < 3; k++)
                weights[f][k] = w[k];
                                    
              wtot += sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
            }

          if(wtot <= 0.0)
            terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);
    
#ifdef SUPERNOVAE     
          double p, Eth;

          if(flag_sn)
            {
              /* Helpers for supernovae injection */
              double num, den, e = 0.0, a = 0.0, b = 0.0;
              double m_ej = WindsAndSN->SN_MassLoss;
              
              /* Ngbs properties (Start with host) */
              int Ngbs = 0;
              double NgbsMass = 0.0, NgbsDensity = 0.0, NgbsMetallicity = 0.0;
            
              double mhost = (P[i].Mass + SphP[i].StarMassFeed);

              Ngbs += 1;
              NgbsMass += mhost * F_HOST;
              NgbsDensity += mhost / SphP[i].Volume * F_HOST;
#ifdef METALS
              NgbsMetallicity += (P[i].Mass * SphP[i].GasMetallicity + SphP[i].StarMetalsFeed) / mhost * F_HOST;
#endif
   
              /* Third pass */ 
              for(f = 0; f < n_faces; f++)
                {
                  q = dc_list[f];

                  int dp = DC[q].dp_index;
                  int particle = Mesh.DP[dp].index;

                  if(particle >= NumGas && Mesh.DP[dp].task == ThisTask)
                    particle -= NumGas;

                  double wbar[3]; 
          
                  for(k = 0; k < 3; k++)
                    wbar[k] = weights[f][k] / wtot;

                  double sq_wbar = (wbar[0]*wbar[0] + wbar[1]*wbar[1] + wbar[2]*wbar[2]);
                  double sqrtsq_wbar = sqrt(sq_wbar);  

                  double mj, vj[3];
                  if(Mesh.DP[dp].task == ThisTask)
                    {
                      mj = P[particle].Mass + SphP[i].StarMassFeed;
                      
                      for(k = 0; k < 3; k++)
                        vj[k] = (SphP[particle].Momentum[k] + SphP[particle].StarMomentumFeed[k]) / mj;   
                    }
                  else
                    {
                      mj = PrimExch[particle].Density * PrimExch[particle].Volume + PrimExch[particle].MassFeed;
                      
                      for(k = 0; k < 3; k++)
                        vj[k] = (PrimExch[particle].Density * PrimExch[particle].Volume * PrimExch[particle].VelGas[k]
                        + PrimExch[particle].MomentumFeed[k]) / mj;

                    }

                  double sq_vj = vj[0]*vj[0] + vj[1]*vj[1] + vj[2]*vj[2];
                  
                  double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0]
                  + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1]
                  + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];
                  
                  double cross = 2.0 * (vj[0]*WindsAndSN->StarVelocity[0]
                  + vj[1]*WindsAndSN->StarVelocity[1]
                  + vj[2]*WindsAndSN->StarVelocity[2]);
         
                  num = 0.5 * mj * m_ej*sqrtsq_wbar*(1-F_HOST) * (sq_vj + sq_vstar - cross);
                  den = mj + m_ej*sqrtsq_wbar*(1-F_HOST);
          
                  e += num / den;
               
                  num = sq_wbar;
                  den = mj + m_ej*sqrtsq_wbar*(1-F_HOST);
          
                  a += num / den;

                  num = mj * ((vj[0] - WindsAndSN->StarVelocity[0]) * wbar[0] 
                  + (vj[1] - WindsAndSN->StarVelocity[1]) * wbar[1]
                  + (vj[2] - WindsAndSN->StarVelocity[2]) * wbar[2]);
                  den = mj + m_ej*sqrtsq_wbar*(1-F_HOST);

                  b += num / den;

                  if(Mesh.DP[dp].task == ThisTask)
                    {
                      Ngbs++;
                      NgbsMass += mj * sqrtsq_wbar * (1-F_HOST);
                      NgbsDensity += mj / SphP[particle].Volume * sqrtsq_wbar * (1-F_HOST);
#ifdef METALS
                      NgbsMetallicity += (P[particle].Mass * SphP[particle].GasMetallicity + SphP[particle].StarMetalsFeed)
                      / mj * sqrtsq_wbar * (1-F_HOST);
#endif
                    }
                  else
                    {
                      Ngbs++;
                      NgbsMass += mj * sqrtsq_wbar * (1-F_HOST);
                      NgbsDensity +=  mj / PrimExch[particle].Volume * sqrtsq_wbar * (1-F_HOST);
#ifdef METALS
                      NgbsMetallicity += (PrimExch[particle].Density * PrimExch[particle].Volume * PrimExch[particle].Scalars[METALS_INDEX] 
                      + PrimExch[particle].MetalsFeed) / mj * sqrtsq_wbar * (1-F_HOST);
#endif
                    }  
                }

              SN_compute(ev, h, e, a, b, NgbsDensity, NgbsMetallicity, &p, &Eth);
            }
#endif
      
          /* Host feedback */
#ifdef WINDS
          if(flag_winds)
            {
              SphP[i].StarMassFeed += WindsAndSN->MassLoss * F_HOST;
              All.StarFeedbackLocal[0] += WindsAndSN->MassLoss * F_HOST;
#ifdef METALS
              SphP[i].StarMetalsFeed += WindsAndSN->MetalsLoss * F_HOST;
              All.StarFeedbackLocal[1] += WindsAndSN->MetalsLoss * F_HOST;
#endif
              for(k = 0; k < 3; k++)
                SphP[i].StarMomentumFeed[k] += WindsAndSN->MassLoss * F_HOST * WindsAndSN->StarVelocity[k];

              double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0] 
              + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1] 
              + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];

              double sq_vwind = WindsAndSN->WindMomentum / WindsAndSN->MassLoss 
              * WindsAndSN->WindMomentum / WindsAndSN->MassLoss;

              SphP[i].StarEnergyFeed += 0.5 * WindsAndSN->MassLoss * F_HOST * (sq_vstar + sq_vwind);
              All.StarFeedbackLocal[2] += 0.5 * WindsAndSN->MassLoss * F_HOST * (sq_vstar + sq_vwind);
            }
#endif 

#ifdef SUPERNOVAE
          if(flag_sn)
            {
              SphP[i].StarMassFeed += WindsAndSN->SN_MassLoss * F_HOST;
              All.StarFeedbackLocal[0] += WindsAndSN->SN_MassLoss * F_HOST;
#ifdef METALS
              SphP[i].StarMetalsFeed += WindsAndSN->SN_MetalsLoss * F_HOST;
              All.StarFeedbackLocal[1] += WindsAndSN->SN_MetalsLoss * F_HOST;
#endif
              for(k = 0; k < 3; k++)
                SphP[i].StarMomentumFeed[k] += WindsAndSN->SN_MassLoss * F_HOST * WindsAndSN->StarVelocity[k];

              SphP[i].StarEnergyFeed += Eth * F_HOST;
              All.StarFeedbackLocal[2] += Eth * F_HOST;
            }
#endif 
 
          /* Fourth pass */  
          for(f = 0; f < n_faces; f++)
            {
              q = dc_list[f];

              double wbar[3]; 
          
              for(k = 0; k < 3; k++)
                wbar[k] = weights[f][k] / wtot;

              double sq_wbar = (wbar[0]*wbar[0] + wbar[1]*wbar[1] + wbar[2]*wbar[2]);
              double sqrtsq_wbar = sqrt(sq_wbar);
      
              struct Feedback_Kick Kick = {0};
              Kick.CellIndex = DC[q].index;
           
              /* Mesh ngbs feedback */
#ifdef WINDS
              if(flag_winds)
                {
                  Kick.DeltaMass = WindsAndSN->MassLoss * sqrtsq_wbar * (1-F_HOST);
#ifdef METALS
                  Kick.DeltaMetals = WindsAndSN->MetalsLoss * sqrtsq_wbar * (1-F_HOST);
#endif    
                  for(k = 0; k < 3; k++)
                    Kick.DeltaP[k] = WindsAndSN->MassLoss * sqrtsq_wbar * (1-F_HOST) 
                    * (WindsAndSN->StarVelocity[k] + WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[k]);
          
                  double sq_vstar = WindsAndSN->StarVelocity[0]*WindsAndSN->StarVelocity[0] 
                  + WindsAndSN->StarVelocity[1]*WindsAndSN->StarVelocity[1] 
                  + WindsAndSN->StarVelocity[2]*WindsAndSN->StarVelocity[2];

                  double sq_vwind = WindsAndSN->WindMomentum / WindsAndSN->MassLoss
                  * WindsAndSN->WindMomentum / WindsAndSN->MassLoss; 

                  double cross = 2.0 * (WindsAndSN->StarVelocity[0] * WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[0] 
                  + WindsAndSN->StarVelocity[1] * WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[1] 
                  + WindsAndSN->StarVelocity[2] * WindsAndSN->WindMomentum / WindsAndSN->MassLoss * wbar[2]);

                  Kick.DeltaE = 0.5 * WindsAndSN->MassLoss * sqrtsq_wbar * (1-F_HOST) 
                  * (sq_vstar + sq_vwind * sq_wbar + cross);
                }   
#endif
 
#ifdef SUPERNOVAE
              if(flag_sn)
                {
                  Kick.SN_DeltaMass = WindsAndSN->SN_MassLoss * sqrtsq_wbar * (1-F_HOST);
#ifdef METALS
                  Kick.SN_DeltaMetals = WindsAndSN->SN_MetalsLoss * sqrtsq_wbar * (1-F_HOST);
#endif 
                  for(k = 0; k < 3; k++)
                    Kick.SN_DeltaP[k] = WindsAndSN->SN_MassLoss * sqrtsq_wbar * (1-F_HOST) * WindsAndSN->StarVelocity[k]
                    + p * wbar[k];

                  Kick.SN_DeltaE = Eth * sqrtsq_wbar * (1-F_HOST);
                }
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
        } //for(int h = 0; h < SphP[i].Host; h++)
      
      /* Go to next host */
      ev += SphP[i].Host;
      /* All stars processed: release host slot */        
      SphP[i].Host = 0;
    } //for(int ev = 0; ev < MechanicalFeedbackEvents.NumEvents;)

  /* MPI exchange of remote kick packets via MPI_Alltoallv */
  int *SendCount = mymalloc("FBSendCount", NTask * sizeof(int));
  int *RecvCount = mymalloc("FBRecvCount", NTask * sizeof(int));
  int *SendDisp = mymalloc("FBSendDisp",  NTask * sizeof(int));
  int *RecvDisp = mymalloc("FBRecvDisp",  NTask * sizeof(int));
 
  memset(SendCount, 0, NTask * sizeof(int));
  for(k = 0; k < n_export; k++)
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
  for(k = 0; k < n_export; k++)
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
  for(k = 0; k < n_recv; k++)
    apply_kick(RecvBuf[k].CellIndex, &RecvBuf[k]);
 
  /* Cleanup in reverse allocation order */
  myfree(RecvBuf);
  myfree(SortedExport);
  myfree(RecvDisp);
  myfree(SendDisp);
  myfree(RecvCount);
  myfree(SendCount);
  myfree(ExportBuf);
  myfree(ExportTask);
 
  TIMER_STOP(CPU_STARS_FEEDBACK);
}