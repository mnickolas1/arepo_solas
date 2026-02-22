#include <math.h>

#include "../main/proto.h"

#include "../stars/star_particle.h"

double StarMassBins[NBINS + 1] = 
{
  /* Region A */
  MMIN, 2.0,

  /* Region B */
  4.0, 6.0, 8.0,

  /* Region C: 8–20, Δm = 0.2 (60 bins) */
  8.2, 8.4, 8.6, 8.8, 9.0,
  9.2, 9.4, 9.6, 9.8, 10.0,
  10.2, 10.4, 10.6, 10.8, 11.0,
  11.2, 11.4, 11.6, 11.8, 12.0,
  12.2, 12.4, 12.6, 12.8, 13.0,
  13.2, 13.4, 13.6, 13.8, 14.0,
  14.2, 14.4, 14.6, 14.8, 15.0,
  15.2, 15.4, 15.6, 15.8, 16.0,
  16.2, 16.4, 16.6, 16.8, 17.0,
  17.2, 17.4, 17.6, 17.8, 18.0,
  18.2, 18.4, 18.6, 18.8, 19.0,
  19.2, 19.4, 19.6, 19.8, 20.0,

  /* Region D: 20–40, Δm = 1 (20 bins) */
  21.0, 22.0, 23.0, 24.0, 25.0,
  26.0, 27.0, 28.0, 29.0, 30.0,
  31.0, 32.0, 33.0, 34.0, 35.0,
  36.0, 37.0, 38.0, 39.0, 40.0,

  /* Region E: 40–80, Δm = 2 (20 bins) */
  42.0, 44.0, 46.0, 48.0, 50.0,
  52.0, 54.0, 56.0, 58.0, 60.0,
  62.0, 64.0, 66.0, 68.0, 70.0,
  72.0, 74.0, 76.0, 78.0, 80.0,

  /* Region F: 80–120, Δm = 4 (10 bins) */
  84.0, 88.0, 92.0, 96.0, 100.0,
  104.0, 108.0, 112.0, 116.0, MMAX
};
  
double StarMeanMassInBins[NBINS];

/* Compute integral with the trapezoid method */
double IntegralTrapezoidal(double a, double b, int N, double (*f)(double))
{
  double h = (b - a) / N;
  double sum = 0.5 * (f(a) + f(b));

  for(int i = 1; i < N; i++)
    sum += f(a + i * h);

  return sum * h;
}
/* Unnormalized Kroupa IMF */
double imf_kroupa(double m) 
{
  if(m < 0.1 || m > 100.0) return 0.0;
  
  if(m < 0.5) return pow(m, -1.3);
  
  return pow(m, -2.3);
}

/* Unnormalized Chabrier IMF */
double imf_chabrier(double m) 
{
  if(m < 0.1 || m > 100.0) return 0.0;

  if(m <= 1.0) 
    {
      double mc = 0.08;
      double sigma = 0.69;
      double logm = log10(m);
      double logmc = log10(mc);
      return (1.0 / m) * exp(-pow((logm - logmc), 2) / (2.0 * sigma*sigma));
    } 
  else return pow(m, -2.3);
}

/* Unnormalized Salpeter IMF */
double imf_salpeter(double m) 
{
  if(m < 0.1 || m > 100.0) return 0.0;
  
  return pow(m, -2.35);
}

/* Select IMF */
double imf(double m) 
{
  switch(All.IMF) 
  {
    case 0: return imf_kroupa(m);
    case 1: return imf_chabrier(m);
    case 2: return imf_salpeter(m);
        
    // fallback
    default: return imf_kroupa(m); 
  }
}

/* Wrapper: m * imf(m) */
double m_times_imf(double m) 
{
  return m * imf(m);
}

void setup_mass_bins(void)
{
  int i;
  double m1, m2, numerator, denominator;
  
  for(i = 0; i < NBINS; i++)
    {
      m1 = StarMassBins[i];
      m2 = StarMassBins[i+1];

      numerator = IntegralTrapezoidal(m1, m2, 100, m_times_imf);
      denominator = IntegralTrapezoidal(m1, m2, 100, imf);

      StarMeanMassInBins[i] = numerator / denominator;
    }
}

/* --- CDF table (built once, reused for all draws) --- */
double cdf_masses[N_CDF_BINS + 1];   /* mass values at each node */
double cdf_values[N_CDF_BINS + 1];   /* cumulative probability at each node */

/* Build a numerical CDF table by integrating the IMF over log-spaced masses */
void build_imf_cdf(void)
{
    double log_mmin = log(MMIN);
    double log_mmax = log(MMAX);
    double dlog = (log_mmax - log_mmin) / N_CDF_BINS;

    /* First pass: fill mass nodes and compute unnormalized cumulative sum.
       We integrate ξ(m) dm = ξ(m) * m * d(ln m), so the integrand in log space
       is imf(m) * m. */
    cdf_masses[0] = MMIN;
    cdf_values[0] = 0.0;

    for(int i = 1; i <= N_CDF_BINS; i++)
      {
        double log_m = log_mmin + i * dlog;
        double m = exp(log_m);
        double log_m_prev = log_mmin + (i-1) * dlog;
        double m_prev = exp(log_m_prev);

        cdf_masses[i] = m;

        /* Trapezoid step in log space: integrand is imf(m)*m */
        double f_left = m_times_imf(m_prev);
        double f_right = m_times_imf(m);
        cdf_values[i] = cdf_values[i-1] + 0.5 * (f_left + f_right) * dlog;
      }

    /* Second pass: normalize so CDF runs from 0 to 1 */
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

    return exp(log(cdf_masses[lo]) + t * (log(cdf_masses[hi]) - log(cdf_masses[lo])));;
}

/* Draw masses for a star particle of total mass M_particle */
void sample_star_particle(double m, int *bins)
{
    /* Zero the bins */
    for(int i = 0; i < NBINS; i++) bins[i] = 0;

    double m_remaining = m;

    while(m_remaining > 0) 
      {
        double u = get_random_number_aux();
        double mstar = sample_imf(u);

        m_remaining -= mstar;
        
        /* Find bin with linear search from bottom */
        int bin = 0;
        while(bin < NBINS - 1 && StarMassBins[bin + 1] < mstar) bin++;
        bins[bin]++;
      }
}