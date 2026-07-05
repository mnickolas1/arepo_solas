#include "../main/allvars.h"
#include "../main/proto.h"

#include "../extern/chealpix.h"


double HealpixDirs[MAX_NUM_RAYS][3];
/* 12*NSIDE^2 */
int NRays; 

/* Reference opacity coefficients [cm^2/g of gas] at solar metallicity */
double Kappa[WAVEBANDS] = 
{
  /* Per gram of GAS at solar Z (= kappa_dust * dust_to_gas_ratio_solar) */
  /* Representative wavelengths: IR~2um, OPT~0.51um, UV~0.15um */
  /* Draine & Li (2007) MW dust model, Cardelli+89 extinction law */
  [INFRARED] = 3.0, /* mid/far-IR dust, cm^2/g, solar Z  */
  [OPTICAL] = 25.0, /* optical-NIR dust, cm^2/g, solar Z */
  [ULTRAVIOLET] = 180.0, /* FUV/photoelectric dust, cm^2/g, solar Z */
  [LYMAN_WERNER] = 0.0, /* DISSOCIATING H2 -> Computed directly from H2 */
  [IONIZING_HI] = 0.0, /* IONIZING_HI -> Computed directly from HI */
  [IONIZING_HeI] = 0.0, /* IONIZING_HeI -> Computed directly from HeI */
  [IONIZING_HeII] = 0.0, /* IONIZING_HeII -> Computed directly from HeII */
};

double ReradiatedFraction[WAVEBANDS] = 
{
  [INFRARED] = 1.0,
  [OPTICAL] = 1.0,
  [ULTRAVIOLET] = 0.95, /* 5% goes to pe heating */
  [LYMAN_WERNER] = 0.0, /* No reradiation->dissociation instead */
  [IONIZING_HI] = 0.0, /* No reradiation->photoionization instead */
  [IONIZING_HeI] = 0.0, /* No reradiation->photoionization instead */
  [IONIZING_HeII] = 0.0, /* No reradiation->photoionization instead */
};

void update_kappa(void)
{
  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;
      
#ifdef METALS
      double Zsol = SphP[i].GasMetallicity / SOLAR_METALLICITY;
#else
      double Zsol = 0;
#endif
      
      double units = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm / All.cf_UnitMass_in_g;
      
      SphP[i].Kappa[INFRARED] = (Kappa[INFRARED] / units) * Zsol;
      SphP[i].Kappa[OPTICAL] = (Kappa[OPTICAL] / units) * Zsol;
      SphP[i].Kappa[ULTRAVIOLET] = (Kappa[ULTRAVIOLET] / units) * Zsol;      

      const double sigma_H2 = 2.47e-18; /* cm^2, effective Solomon-process cross section (Baczynski+15) */ 
      const double sigma_HI = 6.30e-18; /* cm^2, HI threshold, 13.6 eV */
      const double sigma_HeI = 7.83e-18; /* cm^2, HeI threshold, 24.6 eV */
      const double sigma_HeII = 1.58e-18; /* cm^2, HeII threshold, 54.4 eV */
    
      SphP[i].Kappa[LYMAN_WERNER] = (sigma_H2 / (2.0 * PROTONMASS) / units) * SphP[i].GrackleSpecies(GRACKLE_H2I);
      SphP[i].Kappa[IONIZING_HI] = (sigma_HI / PROTONMASS / units) * SphP[i].GrackleSpecies(GRACKLE_HI);
      SphP[i].Kappa[IONIZING_HeI] = (sigma_HeI / (4 * PROTONMASS) / units) * SphP[i].GrackleSpecies(GRACKLE_HeI);
      SphP[i].Kappa[IONIZING_HeII] = (sigma_HeII / (4 * PROTONMASS) / units) * SphP[i].GrackleSpecies(GRACKLE_HeII);  
    }
}

void start_healpix(void) 
{
  int nside = NSIDE_MIN;
  NRays = 12 * nside * nside;

  for(int ipix = 0; ipix < NRays; ipix++)
    {
      pix2vec_nest(nside, ipix, HealpixDirs[ipix]);
    }
}

static RayWorkStack *init_work_stack(long long capacity)
{
  RayWorkStack *w = mymalloc_movable(&w, "RayWorkStack", sizeof(RayWorkStack));
  w->rays = mymalloc_movable(&w->rays, "WorkRays", capacity * sizeof(RayPacket));
  w->n = 0;
  w->capacity = capacity;
  return w;
}

void append_ray(RayWorkStack *w, const RayPacket *ray)
{
  if(w->n >= w->capacity)
    {
      w->capacity *= 2;
      w->rays = myrealloc_movable(w->rays, w->capacity * sizeof(RayPacket));
    }
  w->rays[w->n++] = *ray;
}

static void free_work_stack(RayWorkStack *w)
{
  myfree_movable(w->rays);
  myfree_movable(w);
}

static void init_rays(RayWorkStack *work)
{
  double SQRT3 = sqrt(3);
    
  int ray_idx = 0;
     
  /* Act on host cells */
  for(int ev = 0; ev < MechanicalFeedbackEvents.NumEvents; ev++)
    {
      Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev];
      Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;
        
      /* Loop over rays for this star */
      for(int iray = 0; iray < NRays; iray++)
        {  
          /* Initialize ray from star i */
          RayPacket ray = {0};

          ray.pos[0] = MechanicalFeedback->StarPosition[0];      
          ray.pos[1] = MechanicalFeedback->StarPosition[1]; 
          ray.pos[2] = MechanicalFeedback->StarPosition[2]; 
          ray.dir[0] = HealpixDirs[iray][0];        
          ray.dir[1] = HealpixDirs[iray][1];
          ray.dir[2] = HealpixDirs[iray][2];
          ray.t = 0.0;
          ray.t_exit = MAX_REAL_NUMBER;
          ray.t_maximum = SQRT3 * All.BoxSize;

          ray.active_bands = NO_IR_ACTIVE;

          for(int w = 0; w < WAVEBANDS; w++)
            { 
              ray.Radiated[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
              ray.Radiated[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;

              ray.Radiated_Init[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
              ray.Radiated_Init[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;

              if(ray.Radiated[w].Energy <= 0.0 && ray.Radiated[w].Photons <= 0.0)
                ray.active_bands &= (uint8_t)(~(1u << w));
            }
          
          ray.ray_id = ray_idx;
          ray.home_task = ThisTask;
          
          ray.n_pending = 0;
          ray.target_node = -1;

          ray.is_paused = 0;

          ray.nside = NSIDE_MIN;        
          ray.healpix_pixel = iray;            
            
          ray_idx++;

          append_ray(work, &ray);
        }
    }
}

/* Returns 4 child rays by value - returns 0 if at NSIDE_MAX */
int split_ray(const RayPacket *parent, RayPacket children[4])
{
  if(parent->nside >= NSIDE_MAX)
    return 0;

  int new_nside = parent->nside * 2;

  for(int k = 0; k < 4; k++)
    {
      /* Copy all state including t, active_bands etc */
      children[k] = *parent;   

      children[k].nside = new_nside;
      children[k].healpix_pixel = 4 * parent->healpix_pixel + k;

      pix2vec_nest(new_nside, children[k].healpix_pixel, children[k].dir);

      for(int w = 0; w < WAVEBANDS; w++)
        {
          children[k].Radiated[w].Energy = parent->Radiated[w].Energy * 0.25;
          children[k].Radiated[w].Photons = parent->Radiated[w].Photons * 0.25;

          children[k].Radiated_Init[w].Energy = parent->Radiated_Init[w].Energy * 0.25;
          children[k].Radiated_Init[w].Photons = parent->Radiated_Init[w].Photons * 0.25;
        }
    }
  return 1;
}

static RayExportBuffer *init_export_buffer(int capacity)
{
  RayExportBuffer *buf = malloc(sizeof(RayExportBuffer));
  buf->rays = malloc(capacity * sizeof(RayPacket));
  buf->task = malloc(capacity * sizeof(int));
  buf->n = 0;
  buf->capacity = capacity;
  return buf;
}

static void free_export_buffer(RayExportBuffer *buf)
{
  free(buf->task);
  free(buf->rays);
  free(buf);
}

static void sort_by_task(RayExportBuffer *buf)
{
  /* Build index array */
  int *idx = malloc(buf->n * sizeof(int));
  RayPacket *sorted_rays = malloc(buf->n * sizeof(RayPacket));
  int *sorted_task = malloc(buf->n * sizeof(int));
  
  for(int i = 0; i < buf->n; i++)
    idx[i] = i;

  /* Insertion sort on idx by task */
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

  /* Apply permutation to both arrays */
  for(int i = 0; i < buf->n; i++)
    {
      sorted_rays[i] = buf->rays[idx[i]];
      sorted_task[i] = buf->task[idx[i]];
    }

  memcpy(buf->rays, sorted_rays, buf->n * sizeof(RayPacket));
  memcpy(buf->task, sorted_task, buf->n * sizeof(int));

  free(sorted_task); 
  free(sorted_rays);
  free(idx);
}

static void exchange_rays(RayExportBuffer *send, RayWorkStack *work)
{
  int send_count[NTask], recv_count[NTask];
  int send_offset[NTask], recv_offset[NTask];

  /* Count how many rays go to each task */
  memset(send_count, 0, NTask * sizeof(int));
  for(int i = 0; i < send->n; i++)
    send_count[send->task[i]]++;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  /* Compute offsets */
  send_offset[0] = recv_offset[0] = 0;
  int total_recv = recv_count[0];
    
  for(int i = 1; i < NTask; i++)
    {
      send_offset[i] = send_offset[i-1] + send_count[i-1];
      recv_offset[i] = recv_offset[i-1] + recv_count[i-1];
      total_recv += recv_count[i];
    }

  /* Grow work array to fit incoming rays */
  while(work->n + total_recv > work->capacity)
    {
      work->capacity *= 2;
      work->rays = myrealloc_movable(work->rays, work->capacity * sizeof(RayPacket));
    }

  for(int i = 0; i < NTask; i++) 
    {
      send_count[i] *= sizeof(RayPacket);
      recv_count[i] *= sizeof(RayPacket);
      send_offset[i] *= sizeof(RayPacket);
      recv_offset[i] *= sizeof(RayPacket);
    }

  /* Sort send buffer by task */
  sort_by_task(send);

  /* Exchange ray data */
  MPI_Alltoallv(send->rays, send_count, send_offset, MPI_BYTE, work->rays + work->n, recv_count, recv_offset, MPI_BYTE, MPI_COMM_WORLD);

  work->n += total_recv;
}

#ifdef TREEPOINTS
struct rad_resultsactiveimported_data *Rad_ResultsActiveImported;

static void send_results_home(void)
{
  int i, j, n, k, ncount;
  int *Recv_count = malloc(NTask * sizeof(int));
  int *Send_count = malloc(NTask * sizeof(int));
  int *Recv_offset = malloc(NTask * sizeof(int));
  int *Send_offset = malloc(NTask * sizeof(int));

  /* Count gas cells among imported particles */
  for(i = 0, ncount = 0; i < Tree_NumPartImported; i++)
    if(Tree_Points[i].Type == 0)
      ncount++;

  Rad_ResultsActiveImported = malloc(ncount * sizeof(struct rad_resultsactiveimported_data));

  /* Pack gas cell results */
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
            {
              Rad_ResultsActiveImported[k].Absorbed[w].Energy = Tree_Points[n].Absorbed[w].Energy;
              Rad_ResultsActiveImported[k].Absorbed[w].Photons = Tree_Points[n].Absorbed[w].Photons;
            }

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

  struct rad_resultsactiveimported_data *tmp_results = malloc(Nexport * sizeof(struct rad_resultsactiveimported_data));

  /* Exchange results back to home ranks */
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

  /* Apply results to local particles */
  for(i = 0; i < Nexport; i++)
    {  
      SphP[tmp_results[i].index].StarMomentumFeed[0] += tmp_results[i].StarMomentumFeed[0];
      SphP[tmp_results[i].index].StarMomentumFeed[1] += tmp_results[i].StarMomentumFeed[1];
      SphP[tmp_results[i].index].StarMomentumFeed[2] += tmp_results[i].StarMomentumFeed[2];
      
      for(int w = 0; w < WAVEBANDS; w++)
        {
          SphP[tmp_results[i].index].Absorbed[w].Energy += tmp_results[i].Absorbed[w].Energy;
          SphP[tmp_results[i].index].Absorbed[w].Photons += tmp_results[i].Absorbed[w].Photons;
        }
    }

  /* Free in reverse allocation order */
  free(tmp_results);
  free(Rad_ResultsActiveImported);
  free(Send_offset); free(Recv_offset);
  free(Send_count); free(Recv_count);
}
#endif

#ifdef RAD_OPENING_ANGLE
static void distribute_node_rad(int no)
{
  /* Quick check - skip empty nodes */
  int has_rad = 0;
        
  for(int w = 0; w < WAVEBANDS; w++)
    if(RtNgb_Nodes[no].Absorbed[w].Energy || RtNgb_Nodes[no].Absorbed[w].Photons > 0) 
      { 
        has_rad = 1; 
        break; 
      }
      
  if(!has_rad) 
    return;
  
  double node_tau[WAVEBANDS];
  for(int w = 0; w < WAVEBANDS; w++)
    node_tau[w] = RtNgb_Nodes[no].volume * RtNgb_Nodes[no].density_kappa[w];

  int child = Ngb_Nodes[no].u.d.nextnode;
  while(child != Ngb_Nodes[no].u.d.sibling && child >= 0)
    {
      /* Leaf particle - deposit directly */
      if(child < Ngb_MaxPart)
        {
          if(P[child].Type != 0 || P[child].Mass == 0 || P[child].ID == 0)
            continue;
          
          for(int w = 0; w < WAVEBANDS; w++)
            {
              if(node_tau[w] > 0)
                {
                  double child_tau = SphP[child].Volume * SphP[child].Density * SphP[child].Kappa[w];
                  double frac = child_tau / node_tau[w];
                  
                  SphP[child].Absorbed[w].Energy += frac * RtNgb_Nodes[no].Absorbed[w].Energy;
                  SphP[child].Absorbed[w].Photons += frac * RtNgb_Nodes[no].Absorbed[w].Photons;
                }
            }
          child = Ngb_Nextnode[child];
        }
      /* Internal node - pass fraction down recursively */
      else if(child < Ngb_MaxPart + Ngb_MaxNodes)
        {
          for(int w = 0; w < WAVEBANDS; w++)
            {
              if(node_tau[w] > 0)
                {
                  double child_tau = RtNgb_Nodes[child].volume * RtNgb_Nodes[child].density_kappa[w];
                  double frac = child_tau / node_tau[w];
                  
                  RtNgb_Nodes[child].Absorbed[w].Energy += frac * RtNgb_Nodes[no].Absorbed[w].Energy;
                  RtNgb_Nodes[child].Absorbed[w].Photons += frac * RtNgb_Nodes[no].Absorbed[w].Photons;
                }
            }
          distribute_node_rad(child);
          child = Ngb_Nodes[child].u.d.sibling;
        }
      /* Pseudo-particle - skip, handled by export */
      else
        {    
          child = Ngb_Nextnode[child - Ngb_MaxNodes];
        }
    }

  for(int w = 0; w < WAVEBANDS; w++)
    RtNgb_Nodes[no].Absorbed[w].Energy = RtNgb_Nodes[no].Absorbed[w].Photons = 0.0;
}
#endif

static void radiation_feedback(void)
{
  /* Photoionization and photoelectric heating here -> we do rad pressure inside the tree walk */
  int i, w;
  
  for(i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      double volume = SphP[i].Volume;
      double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval; //do we want the particle dt?

      /* In cgs */
      double V_cgs = volume * (All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm);
      double dt_cgs = dt * All.cf_UnitTime_in_s;

#ifdef PHOTOELECTRIC_HEATING
      double epsilon_pe = 0.05;

      double E_pe = SphP[i].Absorbed[ULTRAVIOLET].Energy * epsilon_pe * All.cf_UnitEnergy_in_cgs; 
      
      /* Volumetric_heating_rate: docs say erg/(s cm^3), straight CGS, no conversion */
      SphP[i].PE_VolHeatingRate +=  E_pe / dt_cgs / V_cgs;
#endif

#ifdef DISSOCIATION
  /* H2 Dissociation */
      double N_abs_H2 = SphP[i].Absorbed[LYMAN_WERNER].Photons;
      
      double n_H2 = SphP[i].GrackleSpecies(GRACKLE_H2I) * SphP[i].Density / (2 * PROTONMASS / All.cf_UnitMass_in_g);
      SphP[i].H2_DissociationRate += n_H2 > 0 ? (N_abs_H2 / dt/All.cf_hubble_a/All.HubbleParam / volume) / n_H2: 0.0;
#endif

#ifdef PHOTOIONIZATION
      double energy_thresh_HI = 13.6 * ELECTRONVOLT_IN_ERGS;
      double energy_thresh_HeI = 24.6 * ELECTRONVOLT_IN_ERGS;
      double energy_thresh_HeII = 54.4 * ELECTRONVOLT_IN_ERGS;

      double E_abs_HI = SphP[i].Absorbed[IONIZING_HI].Energy * All.cf_UnitEnergy_in_cgs;
      double E_abs_HeI = SphP[i].Absorbed[IONIZING_HeI].Energy * All.cf_UnitEnergy_in_cgs;
      double E_abs_HeII = SphP[i].Absorbed[IONIZING_HeII].Energy * All.cf_UnitEnergy_in_cgs;      

      double N_abs_HI = SphP[i].Absorbed[IONIZING_HI].Photons;
      double N_abs_HeI = SphP[i].Absorbed[IONIZING_HeI].Photons;
      double N_abs_HeII = SphP[i].Absorbed[IONIZING_HeII].Photons;

      /* RT_heating_rate: docs say erg/(s cm^3), straight CGS, no conversion */
      double E_threshold_HI = N_abs_HI * energy_thresh_HI; 
      SphP[i].PI_VolHeatingRate += (E_abs_HI - E_threshold_HI) > 0 ? (E_abs_HI - E_threshold_HI) / dt_cgs / V_cgs : 0.0;

      double E_threshold_HeI = N_abs_HeI * energy_thresh_HeI; 
      SphP[i].PI_VolHeatingRate += (E_abs_HeI - E_threshold_HeI) > 0 ? (E_abs_HeI - E_threshold_HeI) / dt_cgs / V_cgs : 0.0;
      
      double E_threshold_HeII = N_abs_HeII * energy_thresh_HeII; 
      SphP[i].PI_VolHeatingRate += (E_abs_HeII - E_threshold_HeII) > 0 ? (E_abs_HeII - E_threshold_HeII) / dt_cgs / V_cgs : 0.0;

      /* RT_ionization_rate:  1 / (time units) */
      double n_HI = SphP[i].GrackleSpecies(GRACKLE_HI) * SphP[i].Density / (PROTONMASS / All.cf_UnitMass_in_g);
      SphP[i].HI_IonizationRate += n_HI > 0 ? (N_abs_HI / dt/All.cf_hubble_a/All.HubbleParam / volume) / n_HI: 0.0;
            
      double n_HeI = SphP[i].GrackleSpecies(GRACKLE_HeI) * SphP[i].Density / (4 * PROTONMASS / All.cf_UnitMass_in_g);
      SphP[i].HeI_IonizationRate += n_HeI > 0 ? (N_abs_HeI / dt/All.cf_hubble_a/All.HubbleParam / volume) / n_HeI: 0.0;

      double n_HeII = SphP[i].GrackleSpecies(GRACKLE_HeII) * SphP[i].Density / (4 * PROTONMASS / All.cf_UnitMass_in_g);
      SphP[i].HeII_IonizationRate += n_HeII > 0 ? (N_abs_HeII / dt/All.cf_hubble_a/All.HubbleParam / volume) / n_HeII: 0.0;
#endif

      for(w = 0; w < WAVEBANDS; w++)
        SphP[i].Absorbed[w].Energy = SphP[i].Absorbed[w].Photons = 0.0;
    }
}

void star_radiation(void)
{
  TIMER_START(CPU_STARS_RADIATION);

  double t0, t1;

#ifdef RAD_OPENING_ANGLE
  /* Zero accumulator on all nodes before treewalk - important for top nodes! */
  for(int no = Ngb_MaxPart; no < Ngb_MaxPart + Ngb_NumNodes; no++)
    {
      for(int w = 0; w < WAVEBANDS; w++)
        {
          RtNgb_Nodes[no].Absorbed[w].Energy = RtNgb_Nodes[no].Absorbed[w].Photons = 0.0;
        }
    }
#endif

  /* Zero accumulator on all leaves before treewalk */
  for(int i = 0; i < NumGas; i++)
    {
      for(int w = 0; w < WAVEBANDS; w++)
        {
          SphP[i].Absorbed[w].Energy = SphP[i].Absorbed[w].Photons = 0.0;
        }
    }
 
  int n_stars = TimeBinsStar.NActiveParticles;
  int n_rays_local = n_stars * NRays;
  long long n_rays_global;

  sumup_large_ints(1, &n_rays_local, &n_rays_global);

  mpi_printf("STAR_RADIATION: Initialize with %d rays\n", n_rays_global);

  RayWorkStack *work = init_work_stack(16 * n_rays_global);
  RayExportBuffer *export_buf = init_export_buffer(5000 * n_rays_global);
  
  init_rays(work);

  long long n_global;
  int iter = 0;
  do
    {
      t0 = second();

      while(work->n > 0)
        {
          RayPacket ray = work->rays[--work->n];
          raytrace_treewalk(&ray, work, export_buf);
        }

      /* Send rays to remote ranks, receive rays from remote ranks */
      exchange_rays(export_buf, work);

      /* Reset export buffer for this round */
      export_buf->n = 0;

      /* Check if anyone still has rays in flight */
      sumup_large_ints(1, &work->n, &n_global);

      t1 = second();

      iter++;

      if(n_global > 0 && iter > 0)
        mpi_printf("STAR_RADIATION: Rad iteration %3d: need to repeat for %12lld rays. (took %g sec)\n", iter, n_global,
        timediff(t0, t1));
      
    } while(n_global > 0);
    
  //send_results_home();

#ifdef RAD_OPENING_ANGLE
  distribute_node_rad(Ngb_MaxPart);
#endif

  radiation_feedback();

  free_export_buffer(export_buf);
  free_work_stack(work);

  TIMER_STOP(CPU_STARS_RADIATION);
}