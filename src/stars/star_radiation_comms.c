#include <limits.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/*
 * Sparse, neighbour-restricted ray exchange
 *
 * Isend/Irecv over the mesh-neighbour graph only (counts)
 * MPI_Iallreduce (termination, posted early)
 * Isend/Irecv over the neighbour graph, zero-count pairs skipped (payload)
 * MPI_Wait (termination)
 */

#define TAG_RAY_COUNT 30201
#define TAG_RAY_DATA 30202

int RayNgbNTask = 0;
int *RayNgbTask = NULL;
int *RayTaskToNgb = NULL;

double RayTMax = 0.0;

static MPI_Datatype MPI_RAYPACKET = MPI_DATATYPE_NULL;

/* All indexed by neighbour slot, allocated once per star_radiation() call */
static int *SendCount = NULL;
static int *RecvCount = NULL;
static int *SendOffset = NULL;
static int *RecvOffset = NULL;
static int *Cursor = NULL;
static MPI_Request *Req = NULL;

/* Scratch for the stable counting sort; persists across rounds */
static int *SortNgb = NULL;
static RayPacket *SortRays = NULL;
static long long SortCapacity = 0;

#ifdef RT_COMM_STATS
double RayCommTimeCounts = 0.0;
double RayCommTimePayload = 0.0;
double RayCommTimeTerm = 0.0;
long long RayCommMsgSent = 0;
long long RayCommRaysSent = 0;
#endif

/*
 * Walk every local cell's Delaunay connection list and flag the ranks that own
 * a face-defining neighbour
 * This is the exact superset of possible export destinations, 
 * so no ray can ever be handed to a rank outside it - and
 * append_export() terminates loudly if one ever is
 */
void ray_comms_init(void)
{
  if(MPI_RAYPACKET == MPI_DATATYPE_NULL)
    {
      MPI_Type_contiguous(sizeof(RayPacket), MPI_BYTE, &MPI_RAYPACKET);
      MPI_Type_commit(&MPI_RAYPACKET);
    }

  char *sflag = malloc(NTask * sizeof(char));
  char *rflag = malloc(NTask * sizeof(char));

  memset(sflag, 0, NTask * sizeof(char));

  for(int i = 0; i < NumGas; i++)
    {
      int q = SphP[i].first_connection;

      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("RAY_COMM: strange connectivity q=%d MaxNvc=%d cell=%d\n", q, MaxNvc, i);

          int dp = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;

          if(particle >= 0)
            {
              const int t = DC[q].task;

              if(t < 0 || t >= NTask)
                terminate("RAY_COMM: DC[%d].task = %d out of range on cell %d\n", q, t, i);

              if(t != ThisTask)
                sflag[t] = 1;
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
    }

  /*
   * Symmetrise: I must be able to RECEIVE from anyone who can send to me
   * AREPO's face connectivity should already be symmetric (the hydro flux
   * exchange relies on it), so in practice rflag == sflag
   */
  MPI_Alltoall(sflag, 1, MPI_CHAR, rflag, 1, MPI_CHAR, MPI_COMM_WORLD);

  RayNgbNTask = 0;
  for(int t = 0; t < NTask; t++)
    if(sflag[t] || rflag[t])
      RayNgbNTask++;

  const int n = RayNgbNTask > 0 ? RayNgbNTask : 1;

  RayNgbTask = malloc(n * sizeof(int));
  RayTaskToNgb = malloc(NTask * sizeof(int));

  for(int t = 0; t < NTask; t++)
    RayTaskToNgb[t] = -1;

  /* Ascending rank order: this is what makes the receive layout match the
     dense Alltoallv layout exactly */
  int k = 0;
  for(int t = 0; t < NTask; t++)
    if(sflag[t] || rflag[t])
      {
        RayTaskToNgb[t] = k;
        RayNgbTask[k++] = t;
      }

  free(rflag);
  free(sflag);

  SendCount = malloc(n * sizeof(int));
  RecvCount = malloc(n * sizeof(int));
  SendOffset = malloc(n * sizeof(int));
  RecvOffset = malloc(n * sizeof(int));
  Cursor = malloc(n * sizeof(int));
  Req = malloc(2 * n * sizeof(MPI_Request));

#ifdef RT_COMM_STATS
  RayCommTimeCounts = RayCommTimePayload = RayCommTimeTerm = 0.0;
  RayCommMsgSent = RayCommRaysSent = 0;
#endif

  int ngb_max, ngb_sum;
  MPI_Allreduce(&RayNgbNTask, &ngb_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&RayNgbNTask, &ngb_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("STAR_RADIATION: RayPacket = %d B, comm neighbours: mean %.1f, max %d (of %d ranks)\n",
             (int)sizeof(RayPacket), (double)ngb_sum / NTask, ngb_max, NTask);
}

void ray_comms_free(void)
{
  free(Req);
  Req = NULL;
  free(Cursor);
  Cursor = NULL;
  free(RecvOffset);
  RecvOffset = NULL;
  free(SendOffset);
  SendOffset = NULL;
  free(RecvCount);
  RecvCount = NULL;
  free(SendCount);
  SendCount = NULL;

  free(RayTaskToNgb);
  RayTaskToNgb = NULL;
  free(RayNgbTask);
  RayNgbTask = NULL;
  RayNgbNTask = 0;

  free(SortRays);
  SortRays = NULL;
  free(SortNgb);
  SortNgb = NULL;
  SortCapacity = 0;

  /* MPI_RAYPACKET is deliberately kept committed across calls */
}

/*
 * Stable counting sort by neighbour slot, using SendOffset[] as the bucket bases
 */
static void sort_by_ngb(RayExportBuffer *buf)
{
  if(buf->n <= 1)
    return;

  if(buf->n > SortCapacity)
    {
      SortCapacity = buf->n;
      SortNgb = realloc(SortNgb, SortCapacity * sizeof(int));
      SortRays = realloc(SortRays, SortCapacity * sizeof(RayPacket));

      if(!SortNgb || !SortRays)
        terminate("RAY_COMM: out of memory growing sort scratch to %lld rays\n", (long long)SortCapacity);
    }

  for(int k = 0; k < RayNgbNTask; k++)
    Cursor[k] = SendOffset[k];

  for(long long i = 0; i < buf->n; i++)
    {
      const int k = buf->ngbs[i];
      const int pos = Cursor[k]++;

      SortNgb[pos] = k;
      SortRays[pos] = buf->rays[i];
    }

  memcpy(buf->ngbs, SortNgb, buf->n * sizeof(int));
  memcpy(buf->rays, SortRays, buf->n * sizeof(RayPacket));
}

void exchange_rays(RayExportBuffer *send, RayWorkStack *work, long long *n_global)
{
  const int nn = RayNgbNTask;

#ifdef RT_COMM_STATS
  double ta, tb;
  ta = second();
#endif

  if(send->n > (long long)INT_MAX)
    terminate("RAY_COMM: %lld rays to export exceeds MPI int count on task %d\n", send->n, ThisTask);

  memset(SendCount, 0, (nn > 0 ? nn : 1) * sizeof(int));

  for(long long i = 0; i < send->n; i++)
    SendCount[send->ngbs[i]]++;

  /* --- counts, over the neighbour graph only --- */
  int nreq = 0;

  for(int k = 0; k < nn; k++)
    MPI_Irecv(&RecvCount[k], 1, MPI_INT, RayNgbTask[k], TAG_RAY_COUNT, MPI_COMM_WORLD, &Req[nreq++]);

  for(int k = 0; k < nn; k++)
    MPI_Isend(&SendCount[k], 1, MPI_INT, RayNgbTask[k], TAG_RAY_COUNT, MPI_COMM_WORLD, &Req[nreq++]);

  MPI_Waitall(nreq, Req, MPI_STATUSES_IGNORE);

  long long total_recv = 0;

  for(int k = 0; k < nn; k++)
    {
      SendOffset[k] = (k == 0) ? 0 : SendOffset[k - 1] + SendCount[k - 1];
      RecvOffset[k] = (k == 0) ? 0 : RecvOffset[k - 1] + RecvCount[k - 1];

      total_recv += RecvCount[k];
    }

#ifdef RT_COMM_STATS
  tb = second();
  RayCommTimeCounts += timediff(ta, tb);
  ta = tb;
#endif

  /*
   * Termination reduce, posted as soon as the post-exchange local count is
   * known and left to complete underneath the payload movement
   */
  long long n_after = work->n + total_recv;

  MPI_Request term_req;
  MPI_Iallreduce(&n_after, n_global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD, &term_req);

  /* Grow BEFORE posting receives: realloc may move work->rays */
  while(work->n + total_recv > work->capacity)
    {
      work->capacity *= 2;
      work->rays = realloc(work->rays, work->capacity * sizeof(RayPacket));

      if(!work->rays)
        terminate("RAY_COMM: out of memory growing work stack to %lld rays\n", work->capacity);
    }

  sort_by_ngb(send);

  /* --- payload, zero-count pairs skipped --- */
  nreq = 0;

  for(int k = 0; k < nn; k++)
    if(RecvCount[k] > 0)
      MPI_Irecv(work->rays + work->n + RecvOffset[k], RecvCount[k], MPI_RAYPACKET,
                RayNgbTask[k], TAG_RAY_DATA, MPI_COMM_WORLD, &Req[nreq++]);

  for(int k = 0; k < nn; k++)
    if(SendCount[k] > 0)
      {
        MPI_Isend(send->rays + SendOffset[k], SendCount[k], MPI_RAYPACKET,
                  RayNgbTask[k], TAG_RAY_DATA, MPI_COMM_WORLD, &Req[nreq++]);
#ifdef RT_COMM_STATS
        RayCommMsgSent++;
        RayCommRaysSent += SendCount[k];
#endif
      }

  MPI_Waitall(nreq, Req, MPI_STATUSES_IGNORE);

  work->n += total_recv;

#ifdef RT_COMM_STATS
  tb = second();
  RayCommTimePayload += timediff(ta, tb);
  ta = tb;
#endif

  MPI_Wait(&term_req, MPI_STATUS_IGNORE);

#ifdef RT_COMM_STATS
  tb = second();
  RayCommTimeTerm += timediff(ta, tb);
#endif
}