#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../extern/chealpix.h"
#include "../stars/star_radiation.h"

double HealpixDirs[MAX_RAYS][3];
int NRays; // 12 * NSIDE^2

/* Opacity coefficients */
double Kappa[WB_COUNT] = {
  1.0,  // LW
  1.0,  // UV
  1.0,  // OP
  1.0   // IR
};

static void sort_by_task(RayExportBuffer *buf)
{
  /* build index array */
  int *idx = malloc(buf->n * sizeof(int));
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
  RayData *sorted_rays = malloc(buf->n * sizeof(RayData));
  int     *sorted_task = malloc(buf->n * sizeof(int));

  for(int i = 0; i < buf->n; i++)
    {
      sorted_rays[i] = buf->rays[idx[i]];
      sorted_task[i] = buf->task[idx[i]];
    }

  memcpy(buf->rays, sorted_rays, buf->n * sizeof(RayData));
  memcpy(buf->task, sorted_task, buf->n * sizeof(int));

  free(sorted_rays);
  free(sorted_task);
  free(idx);
}

void init_healpix_rays(int nside) //run this inside init()
{
    NRays = 12 * nside * nside;
    for(int ipix=0; ipix<NRays; ipix++)
    {
        pix2vec_nest(nside, ipix, HealpixDirs[ipix]);
    }
}

RayData *init_rays_from_stars(int *n_rays_local)
{
  int n_stars = TimeBinsStar.NActiveParticles;
  int total_rays = n_stars * NRays;
    
  // Allocate memory for all rays
  RayData *rays = (RayData *)malloc(total_rays * sizeof(RayData));
    
  if(rays == NULL) 
    {
      fprintf(stderr, "Error: malloc failed for rays\n");
      return NULL;
    }
    
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
          rays[ray_idx].RAD_Ionizing = SP[i].RAD_Ionizing;
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
  RayExportBuffer *buf = malloc(sizeof(RayExportBuffer));
  buf->rays = malloc(capacity * sizeof(RayData));
  buf->task = malloc(capacity * sizeof(int));
  buf->n = 0;
  buf->capacity = capacity;
  return buf;
}

void free_export_buffer(RayExportBuffer *buf)
{
  free(buf->rays);
  free(buf->task);
  free(buf);
}

void radiation(void)
{
    /* 1. initialize rays from active star particles */
    int n_rays_local = 0;
    RayData *rays = init_rays_from_stars(&n_rays_local);

    /* 2. do initial local walk (mode=0) for all rays */
    RayExportBuffer *export_buf = init_export_buffer(MAX_RAYS_EXCHANGE);

    for(int i = 0; i < n_rays_local; i++)
      raytrace_treewalk(&rays[i], 0, -1, export_buf);

    /* 3. iterate until no more exports globally */
    int n_exports_global;
    do
      {
        /* send rays to remote ranks, receive rays from remote ranks */
        RayData *imported_rays;
        int n_imported;
        exchange_rays(export_buf, &imported_rays, &n_imported);

        /* reset export buffer for this round */
        export_buf->n = 0;

        /* walk imported rays in mode=1 */
        for(int i = 0; i < n_imported; i++)
            raytrace_treewalk(&imported_rays[i], 1, 
                              imported_rays[i].target_node, 
                              export_buf);

        myfree(imported_rays);

        /* check if anyone still has rays in flight */
        MPI_Allreduce(&export_buf->n, &n_exports_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      } while(n_exports_global > 0);
    
    send_results_home();

    free_export_buffer(export_buf);
    myfree(rays);
}

void exchange_rays(RayExportBuffer *send, RayData **recv, int *n_recv)
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
       send_count[i] *= sizeof(RayData);
       recv_count[i] *= sizeof(RayData);
       send_offset[i] *= sizeof(RayData);
       recv_offset[i] *= sizeof(RayData);
      }

    /* sort send buffer by task */
    sort_by_task(send);

    *recv = mymalloc("imported_rays", total_recv * sizeof(RayData));
    *n_recv = total_recv;

    /* exchange ray data */
    MPI_Alltoallv(send->rays, send_count, send_offset, MPI_BYTE, *recv, recv_count, recv_offset, MPI_BYTE, MPI_COMM_WORLD);
}

void send_results_home(void)
{
  int i, j, n, k, ncount;
  int *Recv_count = malloc(NTask * sizeof(int));
  int *Send_count = malloc(NTask * sizeof(int));
  int *Recv_offset = malloc(NTask * sizeof(int));
  int *Send_offset = malloc(NTask * sizeof(int));

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
          Rad_ResultsActiveImported[k].RAD_Ionizing = Tree_Points[n].RAD_Ionizing;
          Rad_ResultsActiveImported[k].index        = Tree_Points[n].index;
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
    SphP[tmp_results[i].index].RAD_Ionizing += tmp_results[i].RAD_Ionizing;

  /* free in reverse allocation order */
  myfree(tmp_results);
  myfree(Rad_ResultsActiveImported);
  free(Send_offset); free(Recv_offset);
  free(Send_count);  free(Recv_count);
}