#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"

/* Compute integral with the trapezoid method */
double IntegralTrapezoidal(double a, double b, int N, double (*f)(double))
{
  double h = (b - a) / N;
  double sum = 0.5 * (f(a) + f(b));

  for(int i = 1; i < N; i++)
    sum += f(a + i * h);

  return sum * h;
}

/* Trapezoid in ln m: f(m) dm = f(m) m d(ln m) */
double LogIntegralTrapezoidal(double a, double b, int N, double (*f)(double))
{
  double la = log(a), h = (log(b) - la) / N;
  double sum = 0.5 * (f(a) * a + f(b) * b);

  for(int i = 1; i < N; i++)
    {
      double m = exp(la + i * h);
      sum += f(m) * m;
    }

  return sum * h;
}

/* Unnormalized Kroupa (2001), Maschberger (2013)
   Continuous three-segment power law with breaks at 0.08 and 0.5 Msun */
double imf_kroupa(double m)
{
  if(m < MMIN || m > MMAX)
    return 0.0;

  if(m < 0.08)
    return pow(m, -0.3);

  if(m < 0.5)
    return 0.08 * pow(m, -1.3);

  return 0.04 * pow(m, -2.3);
}

/* Unnormalized Chabrier (2003), Maschberger (2013)
   Continuous between: Lognormal below 1 Msun - Salpeter slope above */
double imf_chabrier(double m)
{
  const double mc = 0.079;
  const double logmc = log10(mc);

  const double sigma = 0.69;

  if(m < MMIN || m > MMAX)
    return 0.0;

  if(m < 1.0)
    {
      double logm = log10(m);
      double d = logm - logmc;
      return (0.158 / m) * exp(-d * d / (2.0 * sigma * sigma));
    }

  return 0.0441 * pow(m, -2.3);
}

/* Unnormalized Salpeter IMF */
double imf_salpeter(double m)
{
  if(m < MMIN || m > MMAX)
    return 0.0;

  return pow(m, -2.35);
}

/* Select IMF */
double imf(double m)
{
  switch(All.IMF)
    {
      case 0:
        return imf_kroupa(m);
      case 1:
        return imf_chabrier(m);
      case 2:
        return imf_salpeter(m);

      /* Fallback */
      default:
        return imf_kroupa(m);
    }
}

/* Wrapper: m * imf(m) */
double m_times_imf(double m)
{
  return m * imf(m);
}

/* --- CDF table (built once, reused for all draws) --- */
double cdf_masses[N_CDF_BINS + 1]; /* mass values at each node */
double cdf_values[N_CDF_BINS + 1]; /* cumulative probability at each node */

/* Build a numerical CDF table by integrating the IMF over log-spaced masses */
void build_imf_cdf(void)
{
  const double log_mmin = log(MMIN), log_mmax = log(MMAX);
  const double dlog = (log_mmax - log_mmin) / N_CDF_BINS;

  cdf_masses[0] = MMIN;
  cdf_values[0] = 0.0;

  double f_prev = m_times_imf(MMIN);

  for(int i = 1; i <= N_CDF_BINS; i++)
    {
      double m = (i == N_CDF_BINS) ? MMAX : exp(log_mmin + i * dlog);
      double f = m_times_imf(m);

      cdf_masses[i] = m;
      cdf_values[i] = cdf_values[i - 1] + 0.5 * (f_prev + f) * dlog;

      f_prev = f;
    }

  double total = cdf_values[N_CDF_BINS];
  for(int i = 0; i <= N_CDF_BINS; i++)
    cdf_values[i] /= total;
}

/* Invert the CDF at a given u in [0,1] using binary search + linear interpolation */
double sample_imf(double u)
{
  /* Binary search for the interval [i, i+1] straddling u */
  int lo = 0, hi = N_CDF_BINS;
  while(hi - lo > 1)
    {
      int mid = (lo + hi) / 2;
      if(cdf_values[mid] <= u)
        lo = mid;
      else
        hi = mid;
    }

  /* Linear interpolation within the interval */
  double cdf_lo = cdf_values[lo];
  double cdf_hi = cdf_values[hi];
  double t = (cdf_hi > cdf_lo) ? (u - cdf_lo) / (cdf_hi - cdf_lo) : 0.0;

  return exp(log(cdf_masses[lo]) + t * (log(cdf_masses[hi]) - log(cdf_masses[lo])));
}

#if defined(STAR_PARTICLES) && STAR_PARTICLES < 2

/* clang-format off */
double StarMassBins[NBINS + 1] =
{
  /* Below LOWEST_MASS_FEEDBACK: a single bin */
  MMIN, 2.0,

  /* 2-8 Msun: no SNe (winds/radiation, AGB) */
  /* 2-4: 20 bins of 0.1 */
  2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3,
  3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4,
  /* 4-6: 10 bins of 0.2 */
  4.2, 4.4, 4.6, 4.8, 5, 5.2, 5.4, 5.6, 5.8, 6,
  /* 6-8: 40 bins of 0.05 */
  6.05, 6.1, 6.15, 6.2, 6.25, 6.3, 6.35, 6.4, 6.45, 6.5,
  6.55, 6.6, 6.65, 6.7, 6.75, 6.8, 6.85, 6.9, 6.95, 7,
  7.05, 7.1, 7.15, 7.2, 7.25, 7.3, 7.35, 7.4, 7.45, 7.5,
  7.55, 7.6, 7.65, 7.7, 7.75, 7.8, 7.85, 7.9, 7.95, 8,

  /* 8-25 Msun: SN progenitors */
  /* 8-12: 100 bins of 0.04 */
  8.04, 8.08, 8.12, 8.16, 8.2, 8.24, 8.28, 8.32, 8.36, 8.4,
  8.44, 8.48, 8.52, 8.56, 8.6, 8.64, 8.68, 8.72, 8.76, 8.8,
  8.84, 8.88, 8.92, 8.96, 9, 9.04, 9.08, 9.12, 9.16, 9.2,
  9.24, 9.28, 9.32, 9.36, 9.4, 9.44, 9.48, 9.52, 9.56, 9.6,
  9.64, 9.68, 9.72, 9.76, 9.8, 9.84, 9.88, 9.92, 9.96, 10,
  10.04, 10.08, 10.12, 10.16, 10.2, 10.24, 10.28, 10.32, 10.36, 10.4,
  10.44, 10.48, 10.52, 10.56, 10.6, 10.64, 10.68, 10.72, 10.76, 10.8,
  10.84, 10.88, 10.92, 10.96, 11, 11.04, 11.08, 11.12, 11.16, 11.2,
  11.24, 11.28, 11.32, 11.36, 11.4, 11.44, 11.48, 11.52, 11.56, 11.6,
  11.64, 11.68, 11.72, 11.76, 11.8, 11.84, 11.88, 11.92, 11.96, 12,
  /* 12-15: 30 bins of 0.1 */
  12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13,
  13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14,
  14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15,
  /* 15-25: 40 bins of 0.25 */
  15.25, 15.5, 15.75, 16, 16.25, 16.5, 16.75, 17, 17.25, 17.5,
  17.75, 18, 18.25, 18.5, 18.75, 19, 19.25, 19.5, 19.75, 20,
  20.25, 20.5, 20.75, 21, 21.25, 21.5, 21.75, 22, 22.25, 22.5,
  22.75, 23, 23.25, 23.5, 23.75, 24, 24.25, 24.5, 24.75, 25,

  /* 25-120 Msun: direct collapse above 24.98 Msun, except SN above 74.83 Msun at most Z <= 0.01 */
  /* 25-50: 10 bins of 2.5 */
  27.5, 30, 32.5, 35, 37.5, 40, 42.5, 45, 47.5, 50,
  /* 50-120: 14 bins of 5 */
  55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
  105, 110, 115, MMAX
};
/* clang-format on */

double StarMeanMassInBins[NBINS];

void setup_mass_bins(void)
{
  int i;
  double m1, m2, numerator, denominator;

  for(i = 0; i < NBINS; i++)
    {
      if(!(StarMassBins[i + 1] > StarMassBins[i]))
        terminate("StarMassBins not strictly increasing at bin %d (%g, %g): check NBINS", i, StarMassBins[i], StarMassBins[i + 1]);
    }

  if(StarMassBins[NBINS] != MMAX)
    terminate("StarMassBins[NBINS] = %g != MMAX = %g: check NBINS", StarMassBins[NBINS], MMAX);

  for(i = 0; i < NBINS; i++)
    {
      m1 = StarMassBins[i];
      m2 = StarMassBins[i + 1];

      numerator = LogIntegralTrapezoidal(m1, m2, 100, m_times_imf);
      denominator = LogIntegralTrapezoidal(m1, m2, 100, imf);

      StarMeanMassInBins[i] = numerator / denominator;
    }
}

static MyStarBins store_bin_counts(int bin, unsigned int n)
{
  if(n <= BIN_COUNTS_MAX)
    return (MyStarBins)n;

#ifdef STAR_FEEDBACK_ACTIVE
  if(StarMeanMassInBins[bin] > LOWEST_MASS_FEEDBACK)
    terminate("Star mass bin %d (%g-%g Msun) holds %u stars, above the STAR_BINS_BITS=%d limit: %u -Raise STAR_BINS_BITS", 
              bin, StarMassBins[bin], StarMassBins[bin + 1], n, (unsigned int)STAR_BINS_BITS, (unsigned int)BIN_COUNTS_MAX);
#endif

  return (MyStarBins)BIN_COUNTS_MAX;
}
#endif

#if defined(STAR_PARTICLES) && STAR_PARTICLES == 0

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *rng;

double norm;
double bin_imf[NBINS];

void setup_imf_integrals(void)
{
  norm = LogIntegralTrapezoidal(MMIN, MMAX, 1000, m_times_imf);

  for(int i = 0; i < NBINS; i++)
    {
      double m1 = StarMassBins[i];
      double m2 = StarMassBins[i + 1];
      bin_imf[i] = LogIntegralTrapezoidal(m1, m2, 100, imf);
    }
}

/* Draw masses for a star particle of total mass M_particle */
void sample_star_particle(double m, MyStarBins *bins)
{
  for(int i = 0; i < NBINS; i++)
    {
      double lambda = m * (bin_imf[i] / norm);
      bins[i] = store_bin_counts(i, (unsigned int) gsl_ran_poisson(rng, lambda));
    }
}
#endif

#if STAR_PARTICLES == 1
/* Draw masses for a star particle of total mass M_particle */
void sample_star_particle(double m, MyStarBins *bins)
{
  /* Zero the bins */
  for(int i = 0; i < NBINS; i++)
    bins[i] = 0;

  double m_sampled = 0.0;

  while(1)
    {
      double u = get_random_number_aux();
      double mstar = sample_imf(u);

      /* Check if adding this star exceeds the target */
      if(m_sampled + mstar > m)
        {
          /* Accept or reject based on which is closer to m */
          double dist_without = m - m_sampled;
          double dist_with = (m_sampled + mstar) - m;

          if(dist_with < dist_without)
            {
              /* Accept: adding the star is closer to m */
              int bin = 0;
              while(bin < NBINS - 1 && StarMassBins[bin + 1] <= mstar)
                bin++;
              bins[bin] = store_bin_counts(bin, (unsigned int)bins[bin] + 1);
            }
          break;
        }

      m_sampled += mstar;

      /* Find bin with linear search from bottom */
      int bin = 0;
      while(bin < NBINS - 1 && StarMassBins[bin + 1] <= mstar)
        bin++;
      bins[bin] = store_bin_counts(bin, (unsigned int)bins[bin] + 1);
    }
}
#endif