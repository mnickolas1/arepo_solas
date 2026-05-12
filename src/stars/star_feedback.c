#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* 
   star_feedback.c

   Loops over gas cells. Any cell with SphP[i].Host > 0 was marked by
   star_density() as hosting a star, and already carries the feedback
   quantities (MassLoss, WindMomentum, SN_EnergyInject, ...) deposited by
   that star in pass 2.

   For each such host cell:
     - scalars (mass, metals, thermal energy) are applied to the host cell
     - momentum (wind + SN) is distributed to face-sharing Voronoi neighbors
       via DC[], with remote neighbors exchanged via MPI_Alltoallv.
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
#endif /* SUPERNOVAE */

/* 
star_feedback() -> main entry point
*/
void star_feedback(void)
{
  CPU_Step[CPU_MISC] += measure_time();

  /* upper bound on export buffer: NumGas host cells * max Voronoi faces */
  int max_export = NumGas * 16;
  struct FBKick *ExportBuf = mymalloc("FBExportBuf",  max_export * sizeof(struct FBKick));
  int *ExportTask = mymalloc("FBExportTask", max_export * sizeof(int));
  int n_export = 0;

#define MAX_FACES 128

  /* 
  Step 1 — loop over gas cells; act on host cells.
  */
  for(int i = 0; i < NumGas; i++)
    {
      if(SphP[i].Host <= 0) continue;   /* not a host cell - skip */

      /* 
      SCALARS - applied directly to the host cell.
      Mass and metals are deposited here; momentum goes to neighbors. 
      */

#ifdef WINDS
      P[i].Mass += SphP[i].MassLoss;
#ifdef METALS
      SphP[i].Metallicity += SphP[i].MetalsLoss;
#endif
#endif /* WINDS */

#ifdef SUPERNOVAE
      P[i].Mass += SphP[i].SN_MassLoss;
#ifdef METALS
      SphP[i].Metallicity += SphP[i].SN_MetalsLoss;
#endif
      double Eth_SN;
      double p0_SN = compute_p0_SN(i, &Eth_SN);
      SphP[i].Energy += Eth_SN;
#endif /* SUPERNOVAE */

      /* 
      FACE-AREA WEIGHTS - pass A over DC[].
      w_f = A_f * max(0, r_hat · n_hat_f)
      where r_hat points from host cell center to face centroid.        
      */
      int dc_list[MAX_FACES];
      double weights[MAX_FACES];
      int n_faces = 0;
      double wtot = 0.0;

      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("star_feedback: strange connectivity q=%d MaxNvc=%d\n", q, MaxNvc);

            if(DC[q].task >= 0 && DC[q].task < NTask)
              {
                int dp = DC[q].dp_index;
                int vf = DC[q].vf_index;
               
                /* outward face normal: vector from cell generator to neighbor
                generator, normalized (this is the Voronoi face normal) */

                double nx = Mesh.DP[dp].x - P[i].Pos[0];
                double ny = Mesh.DP[dp].y - P[i].Pos[1];
                double nz = Mesh.DP[dp].z - P[i].Pos[2];
                double nn = sqrt(nx*nx + ny*ny + nz*nz);

#ifndef REFLECTIVE_X
                if(nx > boxHalf_X)
                  nx -= boxSize_X;
                if(nx < -boxHalf_X)
                  nx += boxSize_X;
#endif /* #ifndef REFLECTIVE_X */

#ifndef REFLECTIVE_Y
                if(ny > boxHalf_Y)
                  ny -= boxSize_Y;
                if(ny < -boxHalf_Y)
                  ny += boxSize_Y;
#endif /* #ifndef REFLECTIVE_Y */

#ifndef REFLECTIVE_Z
                if(nz > boxHalf_Z)
                  nz -= boxSize_Z;
                if(nz < -boxHalf_Z)
                  nz += boxSize_Z;
#endif /* #ifndef REFLECTIVE_Z */

                double w = 0.0;
                if(nn > 0.0 && Mesh.VF[vf].area > 0.0)
                  {
                    nx /= nn;  ny /= nn;  nz /= nn;

                    /* direction from cell center to face centroid */
                    double fcx = Mesh.VF[vf].cx - P[i].Pos[0];
                    double fcy = Mesh.VF[vf].cy - P[i].Pos[1];
                    double fcz = Mesh.VF[vf].cz - P[i].Pos[2];

                    double fr = sqrt(fcx*fcx + fcy*fcy + fcz*fcz);
                    double proj = (fr > 0.0) ? (fcx*nx + fcy*ny + fcz*nz) / fr : 0.0;
                    w = Mesh.VF[vf].area * fmax(0.0, proj);
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
        {
          terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);
        }

      /* 
         MOMENTUM — pass B over DC[]: build and dispatch kick packets. */
      for(int f = 0; f < n_faces; f++)
        {
          q = dc_list[f];

          double wbar = weights[f] / wtot;
          
          if(wbar <= 0.0) 
            terminate("STAR_FEEDBACK: zero weight for host cell %d\n", i);

          /* outward kick direction: host cell center -> face centroid */
          double dx = Mesh.VF[vf].cx - P[i].Pos[0];
          double dy = Mesh.VF[vf].cy - P[i].Pos[1];
          double dz = Mesh.VF[vf].cz - P[i].Pos[2];

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

          double r = sqrt(dx*dx + dy*dy + dz*dz);
          
          if(r <= 0.0) 
            continue;

          struct FBKick kick = {0};
          kick.CellIndex = DC[q].index;

#ifdef WINDS
          kick.DeltaP_wind[0] = wbar * SphP[i].WindMomentum * dx / r;
          kick.DeltaP_wind[1] = wbar * SphP[i].WindMomentum * dy / r;
          kick.DeltaP_wind[2] = wbar * SphP[i].WindMomentum * dz / r;
#endif

#ifdef SUPERNOVAE
          kick.DeltaP_SN[0] = wbar * p0_SN * dx / r;
          kick.DeltaP_SN[1] = wbar * p0_SN * dy / r;
          kick.DeltaP_SN[2] = wbar * p0_SN * dz / r;
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

      /* 
         CLEAR - reset host-cell feedback flags so they don't fire again. */

      SphP[i].Host = 0;

#ifdef WINDS
      SphP[i].MassLoss     = 0.0;
      SphP[i].WindMomentum = 0.0;
#ifdef METALS
      SphP[i].MetalsLoss   = 0.0;
#endif
#endif /* WINDS */

#ifdef SUPERNOVAE
      SphP[i].SN_MassLoss     = 0.0;
      SphP[i].SN_EnergyInject = 0.0;
#ifdef METALS
      SphP[i].SN_MetalsLoss   = 0.0;
#endif
#endif /* SUPERNOVAE */
    }

  /* 
  Step 2 — exchange remote kick packets via MPI_Alltoallv.
  */
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

  /* sort ExportBuf into task-contiguous order via counting sort */
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

  /* 
  Step 3 — apply received kicks to local cells.
  */
  for(int k = 0; k < n_recv; k++)
    apply_kick(RecvBuf[k].CellIndex, &RecvBuf[k]);

  /* cleanup in reverse allocation order */
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