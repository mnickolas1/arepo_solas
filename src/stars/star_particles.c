#include <stdio.h>
#include <math.h>
#include "Integral.h"

/* Unnormalized Kroupa IMF */
double imf_kroupa(double m) 
{
  if (m < 0.1 || m > 100.0) return 0.0;
  if (m < 0.5) return pow(m, -1.3);
  return pow(m, -2.3);
}

/* Unnormalized Chabrier IMF */
double imf_chabrier(double m) 
{
  if (m < 0.1 || m > 100.0) return 0.0;

  if (m <= 1.0) 
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
  if (m < 0.1 || m > 100.0) return 0.0;
  return pow(m, -2.35);
}

/* Select IMF */
double imf(double m) 
{
  switch (All.IMF) 
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

/* Normalization */
double normalization(double M_star, double mmin, double mmax) 
{
  double denom = IntegralTrapezoidal(mmin, mmax, 1000, m_times_imf);
  return M_star / denom;
}

/* Star counts */
double expected_star_count(double imf_norm, double m1, double m2)
{
  return imf_norm * IntegralTrapezoidal(m1, m2, 1000, imf);
}

/* Poisson sampling */
int poisson_sample(double lambda) 
{
  return gsl_ran_poisson(rng, lambda);
}

/* SN events between t and t+dt */
int sne_events(double imf_norm, double mmin, double mmax, double t, double dt)
{
  int N_SN = 0;
  int nbins = 100;
  double dm = (mmax - mmin) / nbins;

  for (int i = 0; i < nbins; i++) 
    {
      double m1 = mmin + i*dm;
      double m2 = m1 + dm;

      double meanN = expected_star_count(imf_norm, m1, m2);
      double m_mid = 0.5*(m1+m2);

      double tau = lifetime(m_mid);
      if (tau >= t && tau < (t+dt)) 
        {
          N_SN += poisson_sample(meanN);
        }
    }
  return N_SN;
}

/* Winds between t and t+dt */
double winds(double imf_norm, double mmin, double mmax, double t, double dt)
{
  double total = 0.0;
  int nbins = 100;
  double dm = (mmax - mmin) / nbins;

  for (int i = 0; i < nbins; i++) 
    {
      double m1 = mmin + i*dm;
      double m2 = m1 + dm;

      double meanN = expected_star_count(imf_norm, m1, m2);
      int Nstars = poisson_sample(meanN);

      double m_mid = 0.5*(m1+m2);
      double mloss = massloss(m_mid, t, dt);

      total += Nstars * mloss;
    }
  return total;
}



