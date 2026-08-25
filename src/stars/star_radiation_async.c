#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"


/*
 * Asynchronous back end
 *
 * Replaces the bulk-synchronous round structure entirely 
 * There is no barrier and no per-round collective: ranks trace whatever they hold, aggregate
 * exports into per-neighbour message buffers, and service MPI between chunks of local work
 * The walk ends when a distributed snapshot proves no rays
 * exist anywhere and none are in flight
 *
 *   1. Send aggregation: One filling buffer per mesh neighbour, 
 *      drawn from a shared slot pool 
 *      A buffer is posted when it fills (RAY_MSG_MAX packets), 
 *      when the rank runs dry, or when RAY_FLUSH_INTERVAL progress
 *      calls have elapsed, so a lone ray never sits in a partial buffer
 *      starving a peer
 *
 *   2. Pre-posted receives: A pool of persistent MPI_ANY_SOURCE Irecvs, each
 *      sized for a full message 
 *      Completed receives are drained straight onto
 *      the work stack and immediately reposted
 *
 *   3. Termination detection: Mattern's four-counter method over a
 *      nonblocking allreduce
 */

#ifndef RT_COMM_SYNC

#define TAG_RAY_DATA 30202

/* Bail out on a genuine deadlock rather than spinning forever */
#define RAY_SLOT_SPIN_MAX 100000000LL

struct RayCommsAsync
{
  RayWorkStack *work; /* imports are pushed here; LIFO gives them priority */

  /* --- send side --- */
  int nslots;
  RayPacket *sendbuf; /* nslots * RAY_MSG_MAX packets, contiguous */
  int *slot_n; /* fill level of each slot */
  MPI_Request *send_req; /* per slot; MPI_REQUEST_NULL when idle */
  int *free_stack;
  int nfree;
  int *fill; /* neighbour slot -> send slot currently filling, or -1 */

  /* --- receive side --- */
  int nrecv;
  RayPacket *recvbuf; /* nrecv * RAY_MSG_MAX packets, contiguous */
  MPI_Request *recv_req;

  /* --- Testsome scratch, sized max(nslots, nrecv) --- */
  int *done_idx;
  MPI_Status *done_stat;

  /* --- Mattern counters (in rays, not messages) --- */
  long long n_sent;
  long long n_recv;
  int activity; /* set on any local progress; cleared when a snapshot posts */

  MPI_Request term_req;
  long long term_in[3];
  long long term_out[3];
  int quiet; /* consecutive clean snapshots; 2 => terminate */
  long long n_snapshots;

  int flush_countdown;

  /* --- stats --- */
  long long msgs_sent;
  long long rays_sent;
  long long rays_traced;
  long long work_hwm;
  long long slot_stalls;
  double trace_time;
};

#define SENDSLOT(c, s) ((c)->sendbuf + (long long)(s) * RAY_MSG_MAX)
#define RECVSLOT(c, r) ((c)->recvbuf + (long long)(r) * RAY_MSG_MAX)

static int drain_recvs(struct RayCommsAsync *c);
static int reclaim_sends(struct RayCommsAsync *c);

/* --------------------------------------------------------------------- */

static void work_reserve(RayWorkStack *w, long long need)
{
  if(w->n + need <= w->capacity)
    return;

  while(w->n + need > w->capacity)
    w->capacity *= 2;

  w->rays = realloc(w->rays, w->capacity * sizeof(RayPacket));

  if(!w->rays)
    terminate("RAY_ASYNC: out of memory growing work stack to %lld rays\n", w->capacity);
}

RayComms *ray_comms_init(RayWorkStack *work)
{
  ray_neighbours_init();

  struct RayCommsAsync *c = calloc(1, sizeof(struct RayCommsAsync));

  c->work = work;

  const int nn = RayNgbNTask > 0 ? RayNgbNTask : 1;

  /* One filling buffer per neighbour, plus spares so a rank exporting hard in
     one direction is not throttled by its own bookkeeping */
  c->nslots = nn + RAY_SEND_SPARE;

  c->sendbuf = malloc((size_t)c->nslots * RAY_MSG_MAX * sizeof(RayPacket));
  c->slot_n = calloc(c->nslots, sizeof(int));
  c->send_req = malloc(c->nslots * sizeof(MPI_Request));
  c->free_stack = malloc(c->nslots * sizeof(int));
  c->fill = malloc(nn * sizeof(int));

  if(!c->sendbuf)
    terminate("RAY_ASYNC: could not allocate %g MB of send buffers\n",
              (double)c->nslots * RAY_MSG_MAX * sizeof(RayPacket) / (1024.0 * 1024.0));

  for(int s = 0; s < c->nslots; s++)
    {
      c->send_req[s] = MPI_REQUEST_NULL;
      c->free_stack[s] = c->nslots - 1 - s; /* pop low indices first */
    }
  c->nfree = c->nslots;

  for(int k = 0; k < nn; k++)
    c->fill[k] = -1;

  /* Enough pre-posted receives that a burst lands directly rather than in the
     unexpected queue, but capped so memory does not scale with NTask */
  c->nrecv = 2 * nn;
  if(c->nrecv > RAY_RECV_SLOTS_MAX)
    c->nrecv = RAY_RECV_SLOTS_MAX;
  if(c->nrecv < RAY_RECV_SLOTS_MIN)
    c->nrecv = RAY_RECV_SLOTS_MIN;

  c->recvbuf = malloc((size_t)c->nrecv * RAY_MSG_MAX * sizeof(RayPacket));
  c->recv_req = malloc(c->nrecv * sizeof(MPI_Request));

  if(!c->recvbuf)
    terminate("RAY_ASYNC: could not allocate %g MB of receive buffers\n",
              (double)c->nrecv * RAY_MSG_MAX * sizeof(RayPacket) / (1024.0 * 1024.0));

  for(int r = 0; r < c->nrecv; r++)
    MPI_Irecv(RECVSLOT(c, r), RAY_MSG_MAX * (int)sizeof(RayPacket), MPI_BYTE,
              MPI_ANY_SOURCE, TAG_RAY_DATA, MPI_COMM_WORLD, &c->recv_req[r]);

  const int nmax = c->nslots > c->nrecv ? c->nslots : c->nrecv;
  c->done_idx = malloc(nmax * sizeof(int));
  c->done_stat = malloc(nmax * sizeof(MPI_Status));

  c->term_req = MPI_REQUEST_NULL;
  c->quiet = 0;

  /* Conservative: forces the first snapshot round to fail, so termination can
     never be declared before any real work has been observed */
  c->activity = 1;

  c->flush_countdown = RAY_FLUSH_INTERVAL;

#ifdef RT_COMM_STATISTICS
  mpi_printf("STAR_RADIATION: async comm, %d send slots + %d recv slots = %.1f MB/rank, %d packets/msg\n",
             c->nslots, c->nrecv,
             (double)(c->nslots + c->nrecv) * RAY_MSG_MAX * sizeof(RayPacket) / (1024.0 * 1024.0),
             RAY_MSG_MAX);
#endif

  return c;
}

/* --------------------------------------------------------------------- */

static void post_send(struct RayCommsAsync *c, int k)
{
  const int s = c->fill[k];

  if(s < 0 || c->slot_n[s] == 0)
    return;

  MPI_Isend(SENDSLOT(c, s), c->slot_n[s] * (int)sizeof(RayPacket), MPI_BYTE,
            RayNgbTask[k], TAG_RAY_DATA, MPI_COMM_WORLD, &c->send_req[s]);

  /* Counted at post time, which is what makes the Mattern snapshot correct:
     a ray on the wire is already in n_sent and not yet in any n_recv */
  c->n_sent += c->slot_n[s];
  c->rays_sent += c->slot_n[s];
  c->msgs_sent++;

  c->fill[k] = -1;
}

static int reclaim_sends(struct RayCommsAsync *c)
{
  int outcount = 0;

  MPI_Testsome(c->nslots, c->send_req, &outcount, c->done_idx, MPI_STATUSES_IGNORE);

  if(outcount == MPI_UNDEFINED || outcount <= 0)
    return 0;

  for(int i = 0; i < outcount; i++)
    {
      const int s = c->done_idx[i];
      c->slot_n[s] = 0;
      c->free_stack[c->nfree++] = s;
    }

  return outcount;
}

static int drain_recvs(struct RayCommsAsync *c)
{
  int outcount = 0;

  MPI_Testsome(c->nrecv, c->recv_req, &outcount, c->done_idx, c->done_stat);

  if(outcount == MPI_UNDEFINED || outcount <= 0)
    return 0;

  long long got = 0;

  for(int i = 0; i < outcount; i++)
    {
      const int r = c->done_idx[i];

      int nbytes = 0;
      MPI_Get_count(&c->done_stat[i], MPI_BYTE, &nbytes);

      if(nbytes % (int)sizeof(RayPacket))
        terminate("RAY_ASYNC: received %d bytes, not a multiple of %d\n", nbytes, (int)sizeof(RayPacket));

      const long long cnt = nbytes / (long long)sizeof(RayPacket);

      if(cnt > 0)
        {
          /* Safe to realloc: nothing holds a pointer into work->rays across a
             call into the comm layer - the driver copies the packet off the
             stack top by value, and split_ray() fills a local array */
          work_reserve(c->work, cnt);

          memcpy(c->work->rays + c->work->n, RECVSLOT(c, r), (size_t)cnt * sizeof(RayPacket));
          c->work->n += cnt;

          c->n_recv += cnt;
          got += cnt;

          if(c->work->n > c->work_hwm)
            c->work_hwm = c->work->n;
        }

      MPI_Irecv(RECVSLOT(c, r), RAY_MSG_MAX * (int)sizeof(RayPacket), MPI_BYTE,
                MPI_ANY_SOURCE, TAG_RAY_DATA, MPI_COMM_WORLD, &c->recv_req[r]);
    }

  if(got > 0)
    c->activity = 1;

  return (int)got;
}

/*
 * Obtain a free send slot
 * Never waits on sends alone: because transfers are
 * rendezvous, a peer's Isend cannot complete until we complete a matching
 * receive, so a rank waiting only on reclaim_sends() could deadlock against a
 * peer doing the same 
 * Draining receives here breaks the cycle - we always
 * have receives posted, and repost them immediately
 */
static int acquire_slot(struct RayCommsAsync *c)
{
  long long spin = 0;

  while(c->nfree == 0)
    {
      reclaim_sends(c);

      if(c->nfree > 0)
        break;

      drain_recvs(c);

      if(++spin > RAY_SLOT_SPIN_MAX)
        terminate("RAY_ASYNC: task %d stalled waiting for a send slot (%d slots, all in flight)\n",
                  ThisTask, c->nslots);
    }

  if(spin > 0)
    c->slot_stalls++;

  return c->free_stack[--c->nfree];
}

void append_export(RayComms *comm, const RayPacket *ray, int task)
{
  struct RayCommsAsync *c = comm;

  const int k = RayTaskToNgb[task];

  if(k < 0)
    terminate("STAR_RADIATION: export to task %d, which is not a mesh neighbour of task %d\n", task, ThisTask);

  if(c->fill[k] < 0)
    c->fill[k] = acquire_slot(c);

  const int s = c->fill[k];

  SENDSLOT(c, s)[c->slot_n[s]++] = *ray;

  if(c->slot_n[s] >= RAY_MSG_MAX)
    post_send(c, k);
}

void ray_comms_flush(RayComms *comm)
{
  struct RayCommsAsync *c = comm;

  for(int k = 0; k < RayNgbNTask; k++)
    if(c->fill[k] >= 0 && c->slot_n[c->fill[k]] > 0)
      post_send(c, k);

  c->flush_countdown = RAY_FLUSH_INTERVAL;
}

void ray_comms_progress(RayComms *comm)
{
  struct RayCommsAsync *c = comm;

  reclaim_sends(c);
  drain_recvs(c);

  /* Bound the latency of a partially filled buffer: without this, a single ray
     headed for a quiet neighbour can sit unsent while this rank grinds through
     unrelated local work, starving that neighbour */
  if(--c->flush_countdown <= 0)
    ray_comms_flush(comm);
}

static int comm_idle(const struct RayCommsAsync *c)
{
  if(c->work->n > 0)
    return 0;

  for(int k = 0; k < RayNgbNTask; k++)
    if(c->fill[k] >= 0 && c->slot_n[c->fill[k]] > 0)
      return 0; /* rays buffered but not yet handed to MPI */

  return 1;
}

/*
 * Mattern's method.
 *
 * A snapshot reduces (sum sent, sum received, any activity since the previous snapshot) 
 * Because the local counters are sampled at different wall-clock
 * times on different ranks, a single clean snapshot is not obviously
 * sufficient: an exact cancellation between one in-flight message and one
 * counter-skew deficit could mask outstanding work
 * Two consecutive clean snapshots close that hole
 *
 * quiet is a pure function of globally reduced values, so every rank computes
 * the same number and all ranks leave the walk on the same snapshot 
 * That matters: these are collectives, and a rank that exited early would hang its
 * peers on the next one 
 */
static int termination_check(struct RayCommsAsync *c)
{
  if(c->term_req == MPI_REQUEST_NULL)
    {
      if(comm_idle(c))
        {
          c->term_in[0] = c->n_sent;
          c->term_in[1] = c->n_recv;
          c->term_in[2] = c->activity;

          c->activity = 0;
          c->n_snapshots++;

          MPI_Iallreduce(c->term_in, c->term_out, 3, MPI_LONG_LONG, MPI_SUM,
                         MPI_COMM_WORLD, &c->term_req);
        }

      return 0;
    }

  int flag = 0;
  MPI_Test(&c->term_req, &flag, MPI_STATUS_IGNORE);

  if(!flag)
    return 0;

  c->term_req = MPI_REQUEST_NULL;

  if(c->term_out[0] == c->term_out[1] && c->term_out[2] == 0)
    c->quiet++;
  else
    c->quiet = 0;

  return c->quiet >= 2;
}

/* --------------------------------------------------------------------- */

void ray_comms_walk(RayWorkStack *work, RayComms *comm)
{
  struct RayCommsAsync *c = comm;

  long long spins = 0;
  int done = 0;

  if(work->n > c->work_hwm)
    c->work_hwm = work->n;

  while(!done)
    {
      int k = 0;

      /* LIFO: imports pushed by drain_recvs() sit on top and are traced first.
         They carry the largest t, so they split hardest and export soonest -
         draining them early keeps neighbours fed and shortens the chain */
      const double t0 = second();

      while(work->n > 0 && k < RAY_PROGRESS_CHUNK)
        {
          RayPacket ray = work->rays[--work->n];
          raytrace_voronoi(&ray, work, comm);
          k++;
        }

      if(k > 0)
        {
          c->activity = 1;
          c->rays_traced += k;
          c->trace_time += timediff(t0, second());
        }

      ray_comms_progress(comm);

      /* About to go idle: partial buffers must go out, or comm_idle() will
         refuse to post a snapshot and termination is never reached */
      if(work->n == 0)
        ray_comms_flush(comm);

      done = termination_check(c);

      //if(++spins >= RAY_ASYNC_SPIN_WARN && (spins % RAY_ASYNC_SPIN_WARN) == 0)
      //  warn("RAY_ASYNC: task %d has spun %lld times (work=%lld sent=%lld recv=%lld snapshots=%lld)\n",
      //       ThisTask, spins, work->n, c->n_sent, c->n_recv, c->n_snapshots);
    }
}

/* --------------------------------------------------------------------- */

void ray_comms_free(RayComms *comm)
{
  struct RayCommsAsync *c = comm;

  /* Termination proved sum(sent) == sum(recv), so every posted send has been
     matched. Waitall is a formality that makes the buffers safe to release */
  MPI_Waitall(c->nslots, c->send_req, MPI_STATUSES_IGNORE);

  /* Tear down the pre-posted receives 
     Any that actually complete here would
     mean a ray arrived after termination was declared, i.e. a bug in the
     detection - so report loudly rather than silently dropping photons */
  for(int r = 0; r < c->nrecv; r++)
    {
      MPI_Cancel(&c->recv_req[r]);

      MPI_Status st;
      MPI_Wait(&c->recv_req[r], &st);

      int cancelled = 0;
      MPI_Test_cancelled(&st, &cancelled);

      if(!cancelled)
        {
          int nbytes = 0;
          MPI_Get_count(&st, MPI_BYTE, &nbytes);
          terminate("RAY_ASYNC: %d bytes of rays arrived after termination on task %d "
                    "(sent=%lld recv=%lld) - termination detection is wrong\n",
                    nbytes, ThisTask, c->n_sent, c->n_recv);
        }
    }

#ifdef RT_COMM_STATISTICS
  {
    long long lsum[5] = {c->msgs_sent, c->rays_sent, c->rays_traced, c->slot_stalls, c->n_snapshots};
    long long gsum[5];
    MPI_Reduce(lsum, gsum, 5, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    long long lmax[2] = {c->work_hwm, c->rays_traced};
    long long gmax[2];
    MPI_Reduce(lmax, gmax, 2, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

    double tmax, tsum;
    MPI_Reduce(&c->trace_time, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&c->trace_time, &tsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    mpi_printf("STAR_RADIATION: async done | %lld msgs, %lld rays exported (%.1f/msg) | "
               "%lld rays traced, max/rank %lld | queue hwm %lld | %lld slot stalls | %lld snapshots\n",
               gsum[0], gsum[1], gsum[0] ? (double)gsum[1] / gsum[0] : 0.0,
               gsum[2], gmax[1], gmax[0], gsum[3], gsum[4]);

    /* Compare against the S_max measured on the sync path: realised speedup
       well below it means the tuning constants, not the algorithm, are wrong */
    mpi_printf("STAR_RADIATION: async trace time max/rank %g s, mean %g s (imbalance %.2f)\n",
               tmax, tsum / NTask, tsum > 0.0 ? tmax * NTask / tsum : 1.0);
  }
#endif

  free(c->done_stat);
  free(c->done_idx);
  free(c->recv_req);
  free(c->recvbuf);
  free(c->fill);
  free(c->free_stack);
  free(c->send_req);
  free(c->slot_n);
  free(c->sendbuf);
  free(c);

  ray_neighbours_free();
}

#endif /* !RT_COMM_SYNC */