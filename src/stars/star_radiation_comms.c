#include <limits.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"


/*
 * Synchronous back end
 */
#ifdef RT_COMM_SYNC

static MPI_Datatype MPI_RAYPACKET = MPI_DATATYPE_NULL;

/* All indexed by neighbour slot, allocated once per star_radiation() call */
static int *SendCount = NULL;
static int *RecvCount = NULL;
static int *SendOffset = NULL;
static int *RecvOffset = NULL;
static int *Cursor = NULL;
static MPI_Request *Req = NULL;

/* Scratch for the stable counting sort */
static int *SortNgbs = NULL;
static RayPacket *SortRays = NULL;
static long long SortCapacity = 0;

static int NRounds = 0;

#ifdef RT_COMM_STATISTICS
static double TimeCounts = 0.0, TimePayload = 0.0, TimeTerm = 0.0;
static double TraceLocalTotal = 0.0, TraceMaxSum = 0.0;
static long long MsgsSent = 0, RaysSent = 0;
#endif

RayComms *ray_comms_init(RayWorkStack *work)
{
  (void)work; /* the synchronous path receives into the work stack directly */

  ray_neighbours_init();

  if(MPI_RAYPACKET == MPI_DATATYPE_NULL)
    {
      MPI_Type_contiguous(sizeof(RayPacket), MPI_BYTE, &MPI_RAYPACKET);
      MPI_Type_commit(&MPI_RAYPACKET);
    }

  const int n = RayNgbNTask > 0 ? RayNgbNTask : 1;

  SendCount = malloc(n * sizeof(int));
  RecvCount = malloc(n * sizeof(int));
  SendOffset = malloc(n * sizeof(int));
  RecvOffset = malloc(n * sizeof(int));
  Cursor = malloc(n * sizeof(int));
  Req = malloc(2 * n * sizeof(MPI_Request));

  NRounds = 0;

#ifdef RT_COMM_STATISTICS
  TimeCounts = TimePayload = TimeTerm = 0.0;
  TraceLocalTotal = TraceMaxSum = 0.0;
  MsgsSent = RaysSent = 0;
#endif

  RayExportBuffer *buf = malloc(sizeof(RayExportBuffer));

  buf->n = 0;
  buf->capacity = 1024;
  buf->ngbs = malloc(buf->capacity * sizeof(int));
  buf->rays = malloc(buf->capacity * sizeof(RayPacket));

  return buf;
}

void append_export(RayComms *comm, const RayPacket *ray, int task)
{
  RayExportBuffer *buf = comm;

  const int k = RayTaskToNgb[task];

  if(k < 0)
    terminate("append_export(): export to task %d, which is not a mesh neighbour of task %d!\n", task, ThisTask);

  if(buf->n >= buf->capacity)
    {
      buf->capacity *= 2;
      buf->ngbs = realloc(buf->ngbs, buf->capacity * sizeof(int));
      buf->rays = realloc(buf->rays, buf->capacity * sizeof(RayPacket));

      if(!buf->ngbs || !buf->rays)
        terminate("append_export(): out of memory growing export buffer to %lld rays!\n", buf->capacity);
    }

  buf->ngbs[buf->n] = k;
  buf->rays[buf->n] = *ray;
  buf->n++;
}

/*
 * Sort rays by destined task
 */
static void sort_by_ngb(RayExportBuffer *buf)
{
  if(buf->n <= 1)
    return;

  if(buf->n > SortCapacity)
    {
      SortCapacity = buf->n;
      SortNgbs = realloc(SortNgbs, SortCapacity * sizeof(int));
      SortRays = realloc(SortRays, SortCapacity * sizeof(RayPacket));

      if(!SortNgbs || !SortRays)
        terminate("sort_by_ngb(): out of memory growing sort scratch to %lld rays!\n", SortCapacity);
    }

  for(int k = 0; k < RayNgbNTask; k++)
    Cursor[k] = SendOffset[k];

  for(long long i = 0; i < buf->n; i++)
    {
      const int k = buf->ngbs[i];
      const int pos = Cursor[k]++;

      SortNgbs[pos] = k;
      SortRays[pos] = buf->rays[i];
    }

  memcpy(buf->ngbs, SortNgbs, buf->n * sizeof(int));
  memcpy(buf->rays, SortRays, buf->n * sizeof(RayPacket));
}

/*
 * Sparse neighbour exchange with the termination reduce posted as soon as the
 * post-exchange local count is known, so it completes underneath the payload
 * movement rather than serialising behind it
 */
static void exchange_rays(RayExportBuffer *send, RayWorkStack *work, long long *n_global)
{
  const int nn = RayNgbNTask;

#ifdef RT_COMM_STATISTICS
  double ta = second(), tb;
#endif

  if(send->n > (long long)INT_MAX)
    terminate("exchange_rays(): %lld rays to export exceeds MPI int count on task %d!\n", send->n, ThisTask);

  memset(SendCount, 0, (nn > 0 ? nn : 1) * sizeof(int));

  for(long long i = 0; i < send->n; i++)
    SendCount[send->ngbs[i]]++;

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

#ifdef RT_COMM_STATISTICS
  tb = second();
  TimeCounts += timediff(ta, tb);
  ta = tb;
#endif

  long long n_after = work->n + total_recv;

  MPI_Request term_req;
  MPI_Iallreduce(&n_after, n_global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD, &term_req);

  /* Grow BEFORE posting receives: realloc may move work->rays */
  while(work->n + total_recv > work->capacity)
    {
      work->capacity *= 2;
      work->rays = realloc(work->rays, work->capacity * sizeof(RayPacket));

      if(!work->rays)
        terminate("exchange_rays(): out of memory growing work stack to %lld rays!\n", work->capacity);
    }

  sort_by_ngb(send);

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

#ifdef RT_COMM_STATISTICS
        MsgsSent++;
        RaysSent += SendCount[k];
#endif
      }

  MPI_Waitall(nreq, Req, MPI_STATUSES_IGNORE);

  work->n += total_recv;

#ifdef RT_COMM_STATISTICS
  tb = second();
  TimePayload += timediff(ta, tb);
  ta = tb;
#endif

  MPI_Wait(&term_req, MPI_STATUS_IGNORE);

#ifdef RT_COMM_STATISTICS
  tb = second();
  TimeTerm += timediff(ta, tb);
#endif
}

void ray_comms_walk(RayWorkStack *work, RayComms *comm)
{
  RayExportBuffer *send = comm;

  long long n_global;
  int iter = 0;

  do
    {
      const double t0 = second();

      /* Do local work first, then exchange */
      while(work->n > 0)
        {
          RayPacket ray = work->rays[--work->n];
          raytrace_voronoi(&ray, work, comm);
        }

#ifdef RT_COMM_STATISTICS
      /* The imbalance diagnostic: sum-of-max is what the barrier charges you,
         max-of-sum is the floor any barrier-free scheme could reach */
      {
        double dt = timediff(t0, second()), dtmax;
        TraceLocalTotal += dt;
        MPI_Allreduce(&dt, &dtmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        TraceMaxSum += dtmax;
      }
#endif

      exchange_rays(send, work, &n_global);

      send->n = 0;

      iter++;

      if(n_global > 0)
        mpi_printf("STAR_RADIATION: Rad iteration %3d: need to repeat for %12lld rays. (took %g sec)\n",
                   iter, n_global, timediff(t0, second()));

      if(iter > 4 * MAXITER)
        terminate("ray_comms_walk(): %lld rays still in flight after %d iterations!\n", n_global, iter);
    }
  while(n_global > 0);

  NRounds = iter;
}

void ray_comms_free(RayComms *comm)
{
  RayExportBuffer *buf = comm;

#ifdef RT_COMM_STATISTICS
  {
    double loc[3] = {TimeCounts, TimePayload, TimeTerm}, mx[3];
    MPI_Reduce(loc, mx, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    long long lmsg[2] = {MsgsSent, RaysSent}, smsg[2];
    MPI_Reduce(lmsg, smsg, 2, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    double worst_rank;
    MPI_Reduce(&TraceLocalTotal, &worst_rank, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    mpi_printf("STAR_RADIATION: sync done | %d rounds | comm max/rank: counts %g s, payload %g s, term %g s "
               "| %lld msgs, %lld rays (%.1f/msg)\n",
               NRounds, mx[0], mx[1], mx[2], smsg[0], smsg[1],
               smsg[0] ? (double)smsg[1] / smsg[0] : 0.0);

    mpi_printf("STAR_RADIATION: imbalance | Smax = %g s, maxS = %g s, S_max = %.2f\n",
               TraceMaxSum, worst_rank, worst_rank > 0.0 ? TraceMaxSum / worst_rank : 1.0);
  }
#endif

  free(buf->rays);
  free(buf->ngbs);
  free(buf);

  free(SortRays);
  SortRays = NULL;
  free(SortNgbs);
  SortNgbs = NULL;
  SortCapacity = 0;

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

  /* MPI_RAYPACKET is deliberately kept committed across calls */

  ray_neighbours_free();
}

#endif /* RT_COMM_SYNC */