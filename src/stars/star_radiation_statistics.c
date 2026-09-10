#include "../main/allvars.h"
#include "../main/proto.h"

#include <stdarg.h>


RTStatistics RTStatisticsLocal;

const char *RayEndNames[RAY_END_CAUSES] =
{
  [RAY_END_TRUNCATE] = "truncate",
  [RAY_END_ESCAPE] = "escape",
  [RAY_END_RELOCATE] = "relocate",
  [RAY_END_STEPCAP] = "stepcap",
  [RAY_END_TMAX] = "tmax",
};

const char *WavebandNames[WAVEBANDS] =
{
  [INFRARED] = "IR", 
  [OPTICAL] = "OP", 
  [ULTRAVIOLET] = "UV",
  [LYMAN_WERNER] = "LW", 
  [IONIZING_HI] = "HI",
  [IONIZING_HeI] = "HeI", 
  [IONIZING_HeII] = "HeII",
};

/* Photon columns are only meaningful for BandTrackPhotons */
void rt_statistics_init(const RayPacket *ray)
{
  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(!(ray->active_bands & (1u << w)))
        continue;

      RTStatisticsLocal.emitted_E[w] += ray->Radiated[w].Energy;

      if((BandTrackPhotons >> w) & 1u)
        RTStatisticsLocal.emitted_N[w] += ray->Radiated[w].Photons;
    }

  RTStatisticsLocal.n_born++;
}

void rt_statistics_reset(void)
{
  memset(&RTStatisticsLocal, 0, sizeof(RTStatisticsLocal));
}

void rt_statistics_drop(const RayPacket *ray, int w)
{
  RTStatisticsLocal.dropped_E[w] += ray->Radiated[w].Energy;

  if((BandTrackPhotons >> w) & 1u)
    RTStatisticsLocal.dropped_N[w] += ray->Radiated[w].Photons;

  RTStatisticsLocal.n_drop[w] += 1.0;
  RTStatisticsLocal.drop_cells[w] += (double)ray->diag_cells;
}

void rt_statistics_abandon(const RayPacket *ray, int cause)
{
  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(!(ray->active_bands & (1u << w)))
        continue;

      RTStatisticsLocal.abandoned_E[w][cause] += ray->Radiated[w].Energy;

      if((BandTrackPhotons >> w) & 1u)
        RTStatisticsLocal.abandoned_N[w][cause] += ray->Radiated[w].Photons;
    }

  RTStatisticsLocal.n_end[cause]++;
  RTStatisticsLocal.end_cells[cause] += (double)ray->diag_cells;
  RTStatisticsLocal.end_t[cause] += ray->t;
}

static void rt_statistics_reduce_d(double *buf, int n)
{
  double *tmp = malloc(n * sizeof(double));

  MPI_Reduce(buf, tmp, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    memcpy(buf, tmp, n * sizeof(double));

  free(tmp);
}

static void rt_statistics_reduce_ll(long long *buf, int n)
{
  long long *tmp = malloc(n * sizeof(long long));

  MPI_Reduce(buf, tmp, n, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    memcpy(buf, tmp, n * sizeof(long long));

  free(tmp);
}

#define BUDGET_LINE_MAX 512

static char *budget_cat(char *p, char *end, const char *fmt, ...)
{
  if(p >= end - 1)
    return p;

  va_list ap;
  va_start(ap, fmt);
  int n = vsnprintf(p, (size_t)(end - p), fmt, ap);
  va_end(ap);

  if(n < 0)
    return p;

  p += n;

  if(p > end - 1)
    p = end - 1;

  return p;
}

void rt_statistics_report(double walltime)
{
  RTStatistics *b = &RTStatisticsLocal;

  rt_statistics_reduce_d(b->emitted_E, WAVEBANDS);
  rt_statistics_reduce_d(b->emitted_N, WAVEBANDS);
  rt_statistics_reduce_d(b->absorbed_E, WAVEBANDS);
  rt_statistics_reduce_d(b->absorbed_N, WAVEBANDS);
  rt_statistics_reduce_d(b->dropped_E, WAVEBANDS);
  rt_statistics_reduce_d(b->dropped_N, WAVEBANDS);
  rt_statistics_reduce_d(&b->abandoned_E[0][0], WAVEBANDS * RAY_END_CAUSES);
  rt_statistics_reduce_d(&b->abandoned_N[0][0], WAVEBANDS * RAY_END_CAUSES);
  rt_statistics_reduce_d(b->drop_cells, WAVEBANDS);
  rt_statistics_reduce_d(b->n_drop, WAVEBANDS);
  rt_statistics_reduce_d(b->end_cells, RAY_END_CAUSES);
  rt_statistics_reduce_d(b->end_t, RAY_END_CAUSES);
  rt_statistics_reduce_ll(b->n_end, RAY_END_CAUSES);

  long long scal[4] = {b->n_born, b->n_split, b->n_crossing, b->n_skipped};
  rt_statistics_reduce_ll(scal, 4);

  if(ThisTask != 0)
    return;

  char line[BUDGET_LINE_MAX], *p, *end;

  double em_tot = 0.0, ab_tot = 0.0, dr_tot = 0.0, lost_tot = 0.0;
  double aban_tot[RAY_END_CAUSES];

  for(int c = 0; c < RAY_END_CAUSES; c++)
    aban_tot[c] = 0.0;

  for(int w = 0; w < WAVEBANDS; w++)
    {
      em_tot += b->emitted_E[w];
      ab_tot += b->absorbed_E[w];
      dr_tot += b->dropped_E[w];

      for(int c = 0; c < RAY_END_CAUSES; c++)
        aban_tot[c] += b->abandoned_E[w][c];
    }

  for(int c = 0; c < RAY_END_CAUSES; c++)
    lost_tot += aban_tot[c];
  lost_tot += dr_tot;

  mpi_printf("\nSTAR_RADIATION: ===== ray statistics =====\n");
  mpi_printf("STAR_RADIATION: %lld rays born, %lld splits, %lld crossings, %lld deposits skipped, walk %g s\n",
             scal[0], scal[1], scal[2], scal[3], walltime);
  mpi_printf("STAR_RADIATION: emitted %.6e (code), deposited %.4f, discarded %.4f\n",
             em_tot, (em_tot > 0.0) ? ab_tot / em_tot : 0.0,
             (em_tot > 0.0) ? lost_tot / em_tot : 0.0);

  /* Per-band energy table, each row normalised by that band's own emission */
  p = line;
  end = line + BUDGET_LINE_MAX;
  p = budget_cat(p, end, "  %-5s %11s %8s %8s", "band", "emitted", "absorb", "drop");
  for(int c = 0; c < RAY_END_CAUSES; c++)
    p = budget_cat(p, end, " %8s", RayEndNames[c]);
  budget_cat(p, end, " %9s", "closure");
  mpi_printf("STAR_RADIATION: energy, as a fraction of each band's emitted budget\n");
  mpi_printf("%s\n", line);

  for(int w = 0; w < WAVEBANDS; w++)
    {
      const double em = b->emitted_E[w];

      p = line;
      end = line + BUDGET_LINE_MAX;
      p = budget_cat(p, end, "  %-5s %11.4e", WavebandNames[w], em);

      if(em <= 0.0)
        {
          budget_cat(p, end, "        -        -");
          mpi_printf("%s\n", line);
          continue;
        }

      double acc = b->absorbed_E[w] + b->dropped_E[w];

      p = budget_cat(p, end, " %8.5f %8.5f", b->absorbed_E[w] / em, b->dropped_E[w] / em);

      for(int c = 0; c < RAY_END_CAUSES; c++)
        {
          acc += b->abandoned_E[w][c];
          p = budget_cat(p, end, " %8.5f", b->abandoned_E[w][c] / em);
        }

      budget_cat(p, end, " %9.1e", (em - acc) / em);
      mpi_printf("%s\n", line);
    }

  /* Photon table, tracked bands only */
  mpi_printf("STAR_RADIATION: photons, tracked bands only\n");
  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(!((BandTrackPhotons >> w) & 1u) || b->emitted_N[w] <= 0.0)
        continue;

      const double em = b->emitted_N[w];
      double acc = b->absorbed_N[w] + b->dropped_N[w];

      p = line;
      end = line + BUDGET_LINE_MAX;
      p = budget_cat(p, end, "  %-5s %11.4e %8.5f %8.5f", WavebandNames[w], em,
                    b->absorbed_N[w] / em, b->dropped_N[w] / em);

      for(int c = 0; c < RAY_END_CAUSES; c++)
        {
          acc += b->abandoned_N[w][c];
          p = budget_cat(p, end, " %8.5f", b->abandoned_N[w][c] / em);
        }

      budget_cat(p, end, " %9.1e", (em - acc) / em);
      mpi_printf("%s\n", line);
    }

  /* Where each band gives up, in cells since the source */
  p = line;
  end = line + BUDGET_LINE_MAX;
  p = budget_cat(p, end, "  mean cells to band drop:");
  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(b->n_drop[w] > 0.0)
        p = budget_cat(p, end, " %s %.1f", WavebandNames[w], b->drop_cells[w] / b->n_drop[w]);
      else
        p = budget_cat(p, end, " %s -", WavebandNames[w]);
    }
  mpi_printf("%s\n", line);

  /* How rays end */
  for(int c = 0; c < RAY_END_CAUSES; c++)
    {
      if(b->n_end[c] == 0)
        continue;

      mpi_printf("  end %-8s %10lld rays, mean %7.1f cells, mean t %10.4e, %7.4f of emitted E\n",
                 RayEndNames[c], b->n_end[c],
                 b->end_cells[c] / (double)b->n_end[c],
                 b->end_t[c] / (double)b->n_end[c],
                 (em_tot > 0.0) ? aban_tot[c] / em_tot : 0.0);
    }
}