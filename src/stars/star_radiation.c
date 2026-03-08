#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../extern/chealpix.h"
#include "../stars/star_radiation.h"

double HealpixDirs[MAX_NUM_RAYS][3];
int NRays; // 12 * NSIDE^2

/* Reference opacity coefficients [cm² g⁻¹ of gas] at solar metallicity.
 *
 * LYMAN_WERNER: FUV dust opacity (~912-2000 Å), dust-to-gas = 0.01 at solar Z.
 *               Ref: Draine (2003), Weingartner & Draine (2001)
 *
 * ULTRAVIOLET:  NUV dust opacity (~2000-4000 Å).
 *               Ref: Draine (2003)
 *
 * OPTICAL:      Visual band (~4000-8000 Å), close to standard V-band extinction.
 *               Ref: Cardelli, Clayton & Mathis (1989)
 *
 * INFRARED:     Modified blackbody at T_DUST_REF = 20 K, beta = 2.
 *               Ref: Semenov et al. (2003), Planck Collaboration XI (2014)
                 
                 Need to check references!
 */
double Kappa[WAVEBANDS] = {
  0,       /* IONIZING_H_PHOTONS     Computed directly from HI              */
  0,       /* IONIZING               Computed directly from HI              */
  1.0e2,   /* LYMAN_WERNER [cm² g⁻¹] — FUV dust at solar Z                  */
  5.0e1,   /* ULTRAVIOLET  [cm² g⁻¹] — NUV dust at solar Z                  */
  1.0e1,   /* OPTICAL      [cm² g⁻¹] — V-band dust at solar Z               */
  1.0e-1,  /* INFRARED     [cm² g⁻¹] — IR dust at T_ref=20K, beta=2, solar Z */ 
};

void update_kappa(void)
{
  for(int i = 0; i < NumGas; i++)
    { 
      double sigma_H = 6.3e-18; // at at Lyman limit
      double Z = SphP[i].GasMetallicity / SOLAR_METALLICITY;
      double units = All.UnitLength_in_cm * All.UnitLength_in_cm / All.UnitMass_in_g;
      
      SphP[i].Kappa[IONIZING_H_PHOTONS] = (sigma_H / PROTONMASS / units) * SphP[i].grHI; 
      SphP[i].Kappa[IONIZING] = (sigma_H / PROTONMASS / units) * SphP[i].grHI;   
      SphP[i].Kappa[LYMAN_WERNER] = (Kappa[LYMAN_WERNER] / units) * Z; 
      SphP[i].Kappa[ULTRAVIOLET] = (Kappa[ULTRAVIOLET] / units) * Z;  
      SphP[i].Kappa[OPTICAL] = (Kappa[OPTICAL] / units) * Z;  
      SphP[i].Kappa[INFRARED] = (Kappa[INFRARED] / units) * Z;  
    }
}

struct rad_resultsactiveimported_data *Rad_ResultsActiveImported;

void init_healpix_rays(void) 
{
  int nside = All.Nside;
  NRays = 12 * nside * nside;

  for(int ipix = 0; ipix < NRays; ipix++)
    {
      pix2vec_nest(All.Nside, ipix, HealpixDirs[ipix]);
    }
}

RayPacket *init_rays_from_stars(int *n_rays_local)
{
  int n_stars = TimeBinsStar.NActiveParticles;
  int total_rays = n_stars * NRays;
    
  // Allocate memory for all rays
  RayPacket *rays = mymalloc("rays", total_rays * sizeof(RayPacket));
    
  *n_rays_local = total_rays;
    
  int ray_idx = 0;
    
  // Loop over all stars
  for(int idx = 0; idx < TimeBinsStar.NActiveParticles; idx++)
    {
      int i = TimeBinsStar.ActiveParticleList[idx];
        
      // Loop over rays for this star
      for(int iray = 0; iray < NRays; iray++)
        {   
          // Initialize ray from star i
          rays[ray_idx].pos[0] = PPS(i).Pos[0];      
          rays[ray_idx].pos[1] = PPS(i).Pos[1]; 
          rays[ray_idx].pos[2] = PPS(i).Pos[2]; 
          rays[ray_idx].dir[0] = HealpixDirs[iray][0];        
          rays[ray_idx].dir[1] = HealpixDirs[iray][1];
          rays[ray_idx].dir[2] = HealpixDirs[iray][2];

          for(int w = 0; w < WAVEBANDS; w++)
            { 
              rays[ray_idx].RAD[w] = SP[i].LUM[w] / NRays;

              rays[ray_idx].RAD_Initial[w] = SP[i].LUM[w] / NRays;
            }
  
          rays[ray_idx].active_bands = ALL_BANDS_ACTIVE;
          
          rays[ray_idx].ray_id = ray_idx;
          rays[ray_idx].home_task = ThisTask;
          
          rays[ray_idx].n_pending = 0;
          rays[ray_idx].target_node = -1;
          rays[ray_idx].t = 0.0;   
            
          ray_idx++;
        }
    }
  return rays;
}

RayExportBuffer *init_export_buffer(int capacity)
{
  RayExportBuffer *buf = mymalloc("export_buf", sizeof(RayExportBuffer));
  buf->rays = mymalloc("export_buf_rays", capacity * sizeof(RayPacket));
  buf->task = mymalloc("export_buf_task", capacity * sizeof(int));
  buf->n = 0;
  buf->capacity = capacity;
  return buf;
}

void free_export_buffer(RayExportBuffer *buf)
{
  myfree(buf->task);
  myfree(buf->rays);
  myfree(buf);
}

static void sort_by_task(RayExportBuffer *buf)
{
  /* build index array */
  int *idx = mymalloc("idx", buf->n * sizeof(int));
  RayPacket *sorted_rays = mymalloc("sorted_rays", buf->n * sizeof(RayPacket));
  int *sorted_task = mymalloc("sorted_task", buf->n * sizeof(int));
  
  for(int i = 0; i < buf->n; i++)
    idx[i] = i;

  /* insertion sort on idx by task */
  for(int i = 1; i < buf->n; i++)
    {
      int key = idx[i];
      int j = i - 1;
      while(j >= 0 && buf->task[idx[j]] > buf->task[key])
        {
          idx[j+1] = idx[j];
          j--;
        }
      idx[j+1] = key;
    }

  /* apply permutation to both arrays */
  for(int i = 0; i < buf->n; i++)
    {
      sorted_rays[i] = buf->rays[idx[i]];
      sorted_task[i] = buf->task[idx[i]];
    }

  memcpy(buf->rays, sorted_rays, buf->n * sizeof(RayPacket));
  memcpy(buf->task, sorted_task, buf->n * sizeof(int));

  myfree(sorted_task); 
  myfree(sorted_rays);
  myfree(idx);
}

void exchange_rays(RayExportBuffer *send, RayPacket **recv, int *n_recv)
{
    int send_count[NTask], recv_count[NTask];
    int send_offset[NTask], recv_offset[NTask];

    /* count how many rays go to each task */
    memset(send_count, 0, NTask * sizeof(int));
    for(int i = 0; i < send->n; i++)
        send_count[send->task[i]]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

    /* compute offsets */
    send_offset[0] = recv_offset[0] = 0;
    int total_recv = recv_count[0];
    
    for(int i = 1; i < NTask; i++)
      {
        send_offset[i] = send_offset[i-1] + send_count[i-1];
        recv_offset[i] = recv_offset[i-1] + recv_count[i-1];
        total_recv += recv_count[i];
      }

    for(int i = 0; i < NTask; i++) 
      {
       send_count[i] *= sizeof(RayPacket);
       recv_count[i] *= sizeof(RayPacket);
       send_offset[i] *= sizeof(RayPacket);
       recv_offset[i] *= sizeof(RayPacket);
      }

    /* sort send buffer by task */
    sort_by_task(send);

    *recv = mymalloc("imported_rays", total_recv * sizeof(RayPacket));
    *n_recv = total_recv;

    /* exchange ray data */
    MPI_Alltoallv(send->rays, send_count, send_offset, MPI_BYTE, *recv, recv_count, recv_offset, MPI_BYTE, MPI_COMM_WORLD);
}

void send_results_home(void)
{
  int i, j, n, k, ncount;
  int *Recv_count = mymalloc("Recv_count", NTask * sizeof(int));
  int *Send_count = mymalloc("Send_count", NTask * sizeof(int));
  int *Recv_offset = mymalloc("Recv_offset", NTask * sizeof(int));
  int *Send_offset = mymalloc("Send_offset", NTask * sizeof(int));

  /* count gas cells among imported particles */
  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
    if(Tree_Points[i].Type == 0)
      ncount++;

  Rad_ResultsActiveImported = mymalloc("Rad_ResultsActiveImported", ncount * sizeof(struct rad_resultsactiveimported_data));

  /* pack gas cell results, count per task */
  for(j = 0; j < NTask; j++)
    Recv_count[j] = 0;

  for(i = 0, n = 0, k = 0; i < NTask; i++)
    for(j = 0; j < Force_Recv_count[i]; j++, n++)
      if(Tree_Points[n].Type == 0)
        {
          Rad_ResultsActiveImported[k].StarMomentumFeed[0] = Tree_Points[n].StarMomentumFeed[0];
          Rad_ResultsActiveImported[k].StarMomentumFeed[1] = Tree_Points[n].StarMomentumFeed[1];
          Rad_ResultsActiveImported[k].StarMomentumFeed[2] = Tree_Points[n].StarMomentumFeed[2];

          for(int w = 0; w < WAVEBANDS; w++)
            Rad_ResultsActiveImported[k].RAD[w] = Tree_Points[n].RAD[w];

          Rad_ResultsActiveImported[k].index = Tree_Points[n].index;
          Recv_count[i]++;
          k++;
        }

  MPI_Alltoall(Recv_count, 1, MPI_INT, Send_count, 1, MPI_INT, MPI_COMM_WORLD);

  int Nexport = 0, Nimport = 0;
  Send_offset[0] = Recv_offset[0] = 0;
  for(j = 0; j < NTask; j++)
    {
      Nexport += Send_count[j];
      Nimport += Recv_count[j];
      if(j > 0)
        {
          Send_offset[j] = Send_offset[j-1] + Send_count[j-1];
          Recv_offset[j] = Recv_offset[j-1] + Recv_count[j-1];
        }
    }

  struct rad_resultsactiveimported_data *tmp_results =
    mymalloc("tmp_results", Nexport * sizeof(struct rad_resultsactiveimported_data));

  /* exchange results back to home ranks */
  for(int ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      int recvTask = ThisTask ^ ngrp;
      if(recvTask < NTask)
        if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
          MPI_Sendrecv(&Rad_ResultsActiveImported[Recv_offset[recvTask]],
                       Recv_count[recvTask] * sizeof(struct rad_resultsactiveimported_data), MPI_BYTE, recvTask, TAG_RAD,
                       &tmp_results[Send_offset[recvTask]],
                       Send_count[recvTask] * sizeof(struct rad_resultsactiveimported_data), MPI_BYTE, recvTask, TAG_RAD,
                       MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

  /* apply results to local particles */
  for(i = 0; i < Nexport; i++)
    {  
      SphP[tmp_results[i].index].StarMomentumFeed[0] += tmp_results[i].StarMomentumFeed[0];
      SphP[tmp_results[i].index].StarMomentumFeed[1] += tmp_results[i].StarMomentumFeed[1];
      SphP[tmp_results[i].index].StarMomentumFeed[2] += tmp_results[i].StarMomentumFeed[2];
      
      for(int w = 0; w < WAVEBANDS; w++)
        SphP[tmp_results[i].index].RAD[w] += tmp_results[i].RAD[w];
    }

  /* free in reverse allocation order */
  myfree(tmp_results);
  myfree(Rad_ResultsActiveImported);
  myfree(Send_offset); myfree(Recv_offset);
  myfree(Send_count);  myfree(Recv_count);
}

void radiation(void)
{
  /* 0. update cell opacities -> maybe we need to do this earlier in the hydro loop */
  update_kappa();

  /* 1. initialize rays from active star particles */
  int n_rays_local = 0;
  RayPacket *rays = init_rays_from_stars(&n_rays_local);

  /* 2. do initial local walk (mode=0) for all rays */
  RayExportBuffer *export_buf = init_export_buffer(MAX_NUM_RAYS_TO_EXCHANGE);

  for(int i = 0; i < n_rays_local; i++)
    raytrace_treewalk(&rays[i], 0, -1, export_buf);

  /* 3. iterate until no more exports globally */
  int n_exports_global;
  do
    {
      /* send rays to remote ranks, receive rays from remote ranks */
      RayPacket *imported_rays;
      int n_imported;
      exchange_rays(export_buf, &imported_rays, &n_imported);

      /* reset export buffer for this round */
      export_buf->n = 0;

      /* walk imported rays in mode=1 */
      for(int i = 0; i < n_imported; i++)
        raytrace_treewalk(&imported_rays[i], 1, imported_rays[i].target_node, export_buf);

      myfree(imported_rays);

      /* check if anyone still has rays in flight */
      MPI_Allreduce(&export_buf->n, &n_exports_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      
    } while(n_exports_global > 0);
    
  send_results_home();

  free_export_buffer(export_buf);
  myfree(rays);
}

void radiation_feedback(void)
{
  /* Photoionization and photoelectric heating here -> we do rad pressure inside the tree walk */
  int i, w;
  
  for(i = 0; i < NumGas; i++)
    {
      double epsilon_pe = 0.05;
      double energy_thresh = 13.6 * ELECTRONVOLT_IN_ERGS;
      double volume = SphP[i].Volume;
      double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;
      //dt *= All.cf_atime / All.cf_time_hubble_a;

      /* in cgs */
      double V_cgs = volume * (All.UnitLength_in_cm * All.UnitLength_in_cm * All.UnitLength_in_cm);
      double dt_cgs = dt * All.UnitTime_in_s;

      double N_abs = SphP[i].RAD[0];
      
      double E_abs = SphP[i].RAD[1] * All.UnitEnergy_in_cgs;
      double E_pe = (SphP[i].RAD[2] + SphP[i].RAD[3]) * epsilon_pe * All.UnitEnergy_in_cgs; // LW + Ultraviolet

      /* RT_ionization_rate */
      double n_HI = SphP[i].grHI * SphP[i].Density / (PROTONMASS / All.UnitMass_in_g);

      SphP[i].HI_IonizationRate = (N_abs / dt / volume) / n_HI;   // 1 / (time units)

      /* RT_heating_rate: docs say erg s⁻¹ cm⁻³, straight CGS, no conversion */
      double E_threshold = N_abs * energy_thresh;                
      SphP[i].PI_VolHeatingRate = ((E_abs - E_threshold) / dt_cgs / V_cgs);

      /* volumetric_heating_rate: docs say erg s⁻¹ cm⁻³, straight CGS, no conversion */
       SphP[i].PE_VolHeatingRate =  E_pe / dt_cgs / V_cgs;

      for(w = 0; w < WAVEBANDS; w++)
        SphP[i].RAD[w] = 0;
    }
}