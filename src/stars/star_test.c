#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef STAR_SB99
/* Table columns: t mdot_w pdot_w mdot_sn snr (L Q) x WAVEBANDS
   in yr, Msun/yr, Msun km/s/yr, Msun/yr, 1/yr, erg/s, 1/s, for SB99_Toma Msun of stars in [LOWEST_MASS_FEEDBACK, MMAX] */
#define SB99_NCOL (5 + 2 * WAVEBANDS)
#define SB99_MAXROWS 100000

static int SB99_N;
static double SB99_Toma;
static double (*SB99)[SB99_NCOL];

void sb99_load(void)
{
  FILE *fd = fopen("sb99_table.txt", "r");
  if(!fd)
    terminate("STAR_SB99: can't open sb99_table.txt in the run directory");

  SB99_N = 0;
  SB99_Toma = 0;
  SB99 = malloc(SB99_MAXROWS * sizeof(*SB99));

  char line[4096];
  while(fgets(line, sizeof(line), fd))
    {
      if(line[0] == '#')
        {
          sscanf(line, "# TOMA %lf", &SB99_Toma);
          continue;
        }

      char *p = line;
      for(int c = 0; c < SB99_NCOL; c++)
        SB99[SB99_N][c] = strtod(p, &p);

      if(++SB99_N == SB99_MAXROWS)
        terminate("STAR_SB99: more than %d rows in sb99_table.txt", SB99_MAXROWS);
    }
  fclose(fd);

  if(SB99_N < 2 || SB99_Toma <= 0)
    terminate("STAR_SB99: sb99_table.txt has %d rows, TOMA = %g", SB99_N, SB99_Toma);

  mpi_printf("STAR_SB99: %d rows, t = %g - %g yr, TOMA = %g Msun\n", SB99_N, SB99[0][0],
             SB99[SB99_N - 1][0], SB99_Toma);
}

/* Rates at age a (linear in t, first row below the grid, zero past it) times dt, scaled to the particle mass:
   the same point sampling as star_feedback_compute()
   State and TimeToSN are kept from the tables */
void sb99_feedback(Star_Feedback *F, double dt, double m, double a)
{
  double r[SB99_NCOL] = {0};

  if(a <= SB99[0][0])
    memcpy(r, SB99[0], sizeof(r));
  else if(a < SB99[SB99_N - 1][0])
    {
      int lo = 0, hi = SB99_N - 1;
      while(hi - lo > 1)
        {
          int mid = (lo + hi) / 2;
          if(SB99[mid][0] <= a)
            lo = mid;
          else
            hi = mid;
        }

      double f = (a - SB99[lo][0]) / (SB99[hi][0] - SB99[lo][0]);
      for(int c = 0; c < SB99_NCOL; c++)
        r[c] = (1 - f) * SB99[lo][c] + f * SB99[hi][c];
    }

  double scale = m / SB99_Toma;

  /* Composition is not compared: fixed solar fractions */
  const double X = 0.7, Y = 0.286, Z = 0.014;

#ifdef WINDS
  F->MassLoss = r[1] * scale * dt;
  F->WindMomentum = r[2] * scale * dt;
#if GRACKLE_CHEMISTRY >= 1
  F->HLoss = X * F->MassLoss;
  F->HeLoss = Y * F->MassLoss;
#endif
#ifdef METALS
  F->MetalsLoss = Z * F->MassLoss;
#endif
#endif

#ifdef SUPERNOVAE
  F->SN_MassLoss = r[3] * scale * dt;
  F->SN_EnergyInject = r[4] * scale * dt * SN_ENERGY;
#if GRACKLE_CHEMISTRY >= 1
  F->SN_HLoss = X * F->SN_MassLoss;
  F->SN_HeLoss = Y * F->SN_MassLoss;
#endif
#ifdef METALS
  F->SN_MetalsLoss = Z * F->SN_MassLoss;
#endif
#endif

#ifdef STAR_RADIATION_ACTIVE
  for(int w = 0; w < WAVEBANDS; w++)
    {
      F->Radiated[w].Energy = r[5 + 2 * w] * scale * dt * SEC_PER_YEAR;
      F->Radiated[w].Photons = r[6 + 2 * w] * scale * dt * SEC_PER_YEAR;
    }
#endif
}
#endif

#ifdef STAR_FB
/* One file per task (the star can migrate between tasks). SFBBIN lines list the sampled mass bins (written at age 0) */
void star_fb(Star_Feedback *F, int i, double dt, double a)
{
  static FILE *fd = NULL;

  if(!fd)
    {
      char fn[MAXLEN_PATH];
      snprintf(fn, MAXLEN_PATH, "%s/sfblog_%04d.txt", All.OutputDir, ThisTask);
      if(!(fd = fopen(fn, "a")))
        terminate("can't open '%s'", fn);
      fprintf(fd, "# SFBBIN ID bin m_lo m_hi m_mean N\n");
      fprintf(fd, "# SFBLOG ID Time age_yr dt_yr Z TimeBinStar HostHydroBin State wind_mass wind_mom SN_mass SN_energy (E N) x 8\n");
    }

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2
  if(a == 0)
    for(int b = 0; b < NBINS; b++)
      if(SP[i].NumOfStarsInBins[b] > 0)
        fprintf(fd, "SFBBIN %llu %d %.8g %.8g %.10g %u\n", (unsigned long long)PPS(i).ID, b, StarMassBins[b], StarMassBins[b + 1],
                StarMeanMassInBins[b], (unsigned int)SP[i].NumOfStarsInBins[b]);
#endif

  double wind[2] = {0}, sn[2] = {0};
#ifdef WINDS
  wind[0] = F->MassLoss;
  wind[1] = F->WindMomentum;
#endif
#ifdef SUPERNOVAE
  sn[0] = F->SN_MassLoss;
  sn[1] = F->SN_EnergyInject;
#endif
#ifdef METALS
  double z = SP[i].Metallicity;
#else
  double z = 0;
#endif

  fprintf(fd, "SFBLOG %llu %.10e %.10e %.10e %.10e %d %d %d %.8e %.8e %.8e %.8e", (unsigned long long)PPS(i).ID, All.Time, a, dt, z,
          SP[i].TimeBinStar, SP[i].HostHydroBin, F->State, wind[0], wind[1], sn[0], sn[1]);
#ifdef STAR_RADIATION_ACTIVE
  for(int w = 0; w < WAVEBANDS; w++)
    fprintf(fd, " %.8e %.8e", F->Radiated[w].Energy, F->Radiated[w].Photons);
#endif
  fprintf(fd, "\n");
  fflush(fd);
}
#endif