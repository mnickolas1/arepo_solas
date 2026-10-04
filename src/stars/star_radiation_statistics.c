#include "../main/allvars.h"
#include "../main/proto.h"

#include <stdarg.h>

RTStatistics RTStatisticsLocal;

const char *RayEndNames[RAY_END_CAUSES] = {
    [RAY_END_TRUNCATE] = "truncate", [RAY_END_ESCAPE] = "escape", [RAY_END_RELOCATE] = "relocate",
    [RAY_END_STEPCAP] = "stepcap",   [RAY_END_TMAX] = "tmax",
};

const char *WavebandNames[WAVEBANDS] = {
    [INFRARED] = "IR",    [OPTICAL] = "OP",     [ULTRAVIOLET] = "UV",   [LYMAN_WERNER] = "LW",
    [IONIZING_HI] = "HI", [IONIZING_H2] = "H2", [IONIZING_HeI] = "HeI", [IONIZING_HeII] = "HeII",
};

/* Photon columns are only meaningful for BandTrackPhotons */
void rt_statistics_init(const RayPacket *ray)
{
  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(!(ray->active_bands & (1u << w)))
        {
          RTStatisticsLocal.untransported_E[w] += ray->Radiated[w].Energy;
          continue;
        }

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
  const double E = ray->Radiated[w].Energy;

  RTStatisticsLocal.dropped_E[w] += E;

  if((BandTrackPhotons >> w) & 1u)
    RTStatisticsLocal.dropped_N[w] += ray->Radiated[w].Photons;

  RTStatisticsLocal.drop_cells[w] += (double)ray->diag_cells;

  RTStatisticsLocal.n_drop[w] += 1.0;
  if(ray->diag_cells == 0)
    RTStatisticsLocal.n_drop_birth[w] += 1.0;
  else if(ray->diag_cells == 1)
    RTStatisticsLocal.n_drop_first[w] += 1.0;
  else
    RTStatisticsLocal.drop_t[w] += ray->t;
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

#define STATISTICS_LINE_MAX 512
/* Column widths, shared by every header and data row in the report */
#define STAT_BAND   5 /* band name */
#define STAT_EM    11 /* emitted, %.4e */
#define STAT_FRAC   8 /* fractions, %.5f */
#define STAT_CLOSE  9 /* closure, %.1e */
#define STAT_COUNT 10 /* event counts, up to ~1e9 */
#define STAT_DEPTH  9 /* mean cells %.1f / mean t %.3e */
#define STAT_CAUSE  9 /* end-cause name */

static char *statistics_cat(char *p, char *end, const char *fmt, ...)
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

/* One centred field of width n, preceded by the column separator */
static char *statistics_mid(char *p, char *end, const char *s, int n)
{
  int len = (int)strlen(s);

  if(len > n)
    len = n;

  const int left = (n - len) / 2;

  return statistics_cat(p, end, " %*s%.*s%*s", left, "", len, s, n - len - left, "");
}

static void statistics_band_header(void)
{
  char line[STATISTICS_LINE_MAX];
  char *p = line, *end = line + STATISTICS_LINE_MAX;

  p = statistics_cat(p, end, " ");
  p = statistics_mid(p, end, "band", STAT_BAND);
  p = statistics_mid(p, end, "emitted", STAT_EM);
  p = statistics_mid(p, end, "absorb", STAT_FRAC);
  p = statistics_mid(p, end, "drop", STAT_FRAC);

  for(int c = 0; c < RAY_END_CAUSES; c++)
    p = statistics_mid(p, end, RayEndNames[c], STAT_FRAC);

  statistics_mid(p, end, "closure", STAT_CLOSE);

  mpi_printf("%s\n", line);
}

void rt_statistics_report(void)
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
  rt_statistics_reduce_d(b->drop_t, WAVEBANDS);
  rt_statistics_reduce_d(b->n_drop, WAVEBANDS);
  rt_statistics_reduce_d(b->n_drop_birth, WAVEBANDS);
  rt_statistics_reduce_d(b->n_drop_first, WAVEBANDS);
  rt_statistics_reduce_d(b->untransported_E, WAVEBANDS);
  rt_statistics_reduce_d(b->end_cells, RAY_END_CAUSES);
  rt_statistics_reduce_d(b->end_t, RAY_END_CAUSES);
  rt_statistics_reduce_ll(b->n_end, RAY_END_CAUSES);

  long long scal[4] = {b->n_born, b->n_split, b->n_crossing, b->n_skipped};
  rt_statistics_reduce_ll(scal, 4);

  if(ThisTask != 0)
    return;

  char line[STATISTICS_LINE_MAX], *p, *end;

  double em_tot = 0.0, ab_tot = 0.0, dr_tot = 0.0, lo_tot = 0.0, un_tot = 0.0;
  double aband_tot[RAY_END_CAUSES];

  for(int c = 0; c < RAY_END_CAUSES; c++)
    aband_tot[c] = 0.0;

  for(int w = 0; w < WAVEBANDS; w++)
    {
      em_tot += b->emitted_E[w];
      ab_tot += b->absorbed_E[w];
      dr_tot += b->dropped_E[w];
      un_tot += b->untransported_E[w];

      for(int c = 0; c < RAY_END_CAUSES; c++)
        aband_tot[c] += b->abandoned_E[w][c];
    }

  for(int c = 0; c < RAY_END_CAUSES; c++)
    lo_tot += aband_tot[c];
  lo_tot += dr_tot;

  mpi_printf("\nSTAR_RADIATION: ===== ray statistics =====\n");
  mpi_printf("STAR_RADIATION: %lld rays born, %lld splits, %lld crossings, %lld deposits skipped\n",
             scal[0], scal[1], scal[2], scal[3]);
  mpi_printf("STAR_RADIATION: emitted %.6e (code), deposited %.4f, discarded %.4f\n",
             em_tot, (em_tot > 0.0) ? ab_tot / em_tot : 0.0,
             (em_tot > 0.0) ? lo_tot / em_tot : 0.0);

  /* Per-band energy table, each row normalised by that band's own emission */
  mpi_printf("STAR_RADIATION: energy, as a fraction of each band's emitted statistics\n");
  statistics_band_header();

  for(int w = 0; w < WAVEBANDS; w++)
    {
      const double em = b->emitted_E[w];

      p = line;
      end = line + STATISTICS_LINE_MAX;
      p = statistics_cat(p, end, " ");
      p = statistics_mid(p, end, WavebandNames[w], STAT_BAND);
      p = statistics_cat(p, end, " %*.4e", STAT_EM, em);

      if(em <= 0.0)
        {
          for(int k = 0; k < 2 + RAY_END_CAUSES; k++)
            p = statistics_cat(p, end, " %*s", STAT_FRAC, "-");

          statistics_cat(p, end, " %*s", STAT_CLOSE, "-");
          mpi_printf("%s\n", line);
          continue;
        }

      double acc = b->absorbed_E[w] + b->dropped_E[w];

      p = statistics_cat(p, end, " %*.5f %*.5f", STAT_FRAC, b->absorbed_E[w] / em,
                         STAT_FRAC, b->dropped_E[w] / em);

      for(int c = 0; c < RAY_END_CAUSES; c++)
        {
          acc += b->abandoned_E[w][c];
          p = statistics_cat(p, end, " %*.5f", STAT_FRAC, b->abandoned_E[w][c] / em);
        }

      statistics_cat(p, end, " %*.1e", STAT_CLOSE, (em - acc) / em);
      mpi_printf("%s\n", line);
    }

  /* Photon table, tracked bands only */
  mpi_printf("STAR_RADIATION: photons, tracked bands only\n");
  statistics_band_header();

  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(!((BandTrackPhotons >> w) & 1u) || b->emitted_N[w] <= 0.0)
        continue;

      const double em = b->emitted_N[w];
      double acc = b->absorbed_N[w] + b->dropped_N[w];

      p = line;
      end = line + STATISTICS_LINE_MAX;
      p = statistics_cat(p, end, " ");
      p = statistics_mid(p, end, WavebandNames[w], STAT_BAND);
      p = statistics_cat(p, end, " %*.4e %*.5f %*.5f", STAT_EM, em,
                         STAT_FRAC, b->absorbed_N[w] / em,
                         STAT_FRAC, b->dropped_N[w] / em);

      for(int c = 0; c < RAY_END_CAUSES; c++)
        {
          acc += b->abandoned_N[w][c];
          p = statistics_cat(p, end, " %*.5f", STAT_FRAC, b->abandoned_N[w][c] / em);
        }

      statistics_cat(p, end, " %*.1e", STAT_CLOSE, (em - acc) / em);
      mpi_printf("%s\n", line);
    }

  mpi_printf("STAR_RADIATION: band drop depth (* = genuine drops, excl. birth and 1cell)\n");

  p = statistics_cat(p, end, " ");
  p = statistics_mid(p, end, "band", STAT_BAND);
  p = statistics_mid(p, end, "n_drop", STAT_COUNT);
  p = statistics_mid(p, end, "birth", STAT_COUNT);
  p = statistics_mid(p, end, "1cell", STAT_COUNT);
  p = statistics_mid(p, end, "<cells>", STAT_DEPTH);
  p = statistics_mid(p, end, "<cells>*", STAT_DEPTH);
  statistics_mid(p, end, "<t>*", STAT_DEPTH);
  mpi_printf("%s\n", line);

  for(int w = 0; w < WAVEBANDS; w++)
    {
      const double n = b->n_drop[w];
      const double n_deep = n - b->n_drop_birth[w] - b->n_drop_first[w];

      p = line;
      end = line + STATISTICS_LINE_MAX;
      p = statistics_cat(p, end, " ");
      p = statistics_mid(p, end, WavebandNames[w], STAT_BAND);

      if(n <= 0.0)
        {
          for(int k = 0; k < 3; k++)
            p = statistics_cat(p, end, " %*s", STAT_COUNT, "-");

          for(int k = 0; k < 3; k++)
            p = statistics_cat(p, end, " %*s", STAT_DEPTH, "-");

          mpi_printf("%s\n", line);
          continue;
        }

      p = statistics_cat(p, end, " %*.0f %*.0f %*.0f", STAT_COUNT, n,
                         STAT_COUNT, b->n_drop_birth[w],
                         STAT_COUNT, b->n_drop_first[w]);

      p = statistics_cat(p, end, " %*.1f", STAT_DEPTH, b->drop_cells[w] / n);

      if(n_deep > 0.0)
        {
          p = statistics_cat(p, end, " %*.1f", STAT_DEPTH,
                             (b->drop_cells[w] - b->n_drop_first[w]) / n_deep);
          statistics_cat(p, end, " %*.3e", STAT_DEPTH, b->drop_t[w] / n_deep);
        }
      else
        {
          p = statistics_cat(p, end, " %*s", STAT_DEPTH, "-");
          statistics_cat(p, end, " %*s", STAT_DEPTH, "-");
        }

      mpi_printf("%s\n", line);
    }

  if(un_tot > 0.0)
    {
      p = line;
      end = line + STATISTICS_LINE_MAX;
      p = statistics_cat(p, end, "STAR_RADIATION: emitted but never transported:");

      for(int w = 0; w < WAVEBANDS; w++)
        {
          if(b->untransported_E[w] <= 0.0)
            continue;

          p = statistics_cat(p, end, " %s %.4e", WavebandNames[w], b->untransported_E[w]);
        }

      statistics_cat(p, end, " (code), %.1f%% of source output",
                     100.0 * un_tot / (em_tot + un_tot));
      mpi_printf("%s\n", line);
    }

  /* How rays end */
  mpi_printf("STAR_RADIATION: how rays end\n");

  p = line;
  end = line + STATISTICS_LINE_MAX;
  p = statistics_cat(p, end, " ");
  p = statistics_mid(p, end, "cause", STAT_CAUSE);
  p = statistics_mid(p, end, "rays", STAT_COUNT);
  p = statistics_mid(p, end, "<cells>", STAT_DEPTH);
  p = statistics_mid(p, end, "<t>", STAT_DEPTH);
  statistics_mid(p, end, "frac_E", STAT_FRAC);
  mpi_printf("%s\n", line);

  for(int c = 0; c < RAY_END_CAUSES; c++)
    {
      if(b->n_end[c] == 0)
        continue;

      p = line;
      end = line + STATISTICS_LINE_MAX;
      p = statistics_cat(p, end, " ");
      p = statistics_mid(p, end, RayEndNames[c], STAT_CAUSE);
      statistics_cat(p, end, " %*lld %*.1f %*.3e %*.4f", STAT_COUNT, b->n_end[c],
                     STAT_DEPTH, b->end_cells[c] / (double)b->n_end[c],
                     STAT_DEPTH, b->end_t[c] / (double)b->n_end[c],
                     STAT_FRAC, (em_tot > 0.0) ? aband_tot[c] / em_tot : 0.0);
      mpi_printf("%s\n", line);
    }
}