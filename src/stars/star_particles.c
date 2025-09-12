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

void star_particle_feedback()
{  
   double mmin = 8.0;      // min stellar mass
   double mmax = 100.0;    // max stellar mass

   for (int i=0; i<NumStars; i++)
   {
      // 1) Compute IMF normalization for this particle
      double imf_norm = normalization(M_star, mmin, mmax);

      // 2) Count SN events this timestep
      int N_SN = count_sne(imf_norm, mmin, mmax, t, dt);

      // 3) Compute total wind mass loss this timestep
      double mass_loss = winds(imf_norm, mmin, mmax, t, dt);
   }
}

/*  
*  Feedback tables interpolation
*/

/* Linear interpolation helper function */
double linear_interpolation(double x, double x0, double x1, double y0, double y1) 
{
  // avoid divide by zero
  if (x1 == x0) return y0;

  return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

/* Linear interpolation in age */
double interpolate_age(int track, double t) 
{
  const double *ages  = age_arrays[track];
  const double *mloss = mloss_arrays[track];
  int N = nsteps[track];

  if (t <= ages[0])   return mloss[0];
  if (t >= ages[N-1]) return mloss[N-1];

  for (int i = 0; i < N-1; i++) 
  {
    if (t >= ages[i] && t <= ages[i+1]) 
    {
      return linear_interpolation(t, ages[i], ages[i+1], mloss[i], mloss[i+1]);
    }
  }
  // fallback
  return mloss[N-1]; 
}

/* Linear interpolation in stellar mass */
double interpolate_stellar_mass(double Mstar_init, double age) 
{
  //units (years - solar masses)
  age *= All.UnitTime_in_s / SEC_PER_YEAR;
  Mstar_init *= All.UnitMass_in_g / SOLAR_MASS;

  //only include winds of massive stars
  if(Mstar_init < 8) return 0.0;

  if (Mstar_init <= init_mass[0])
    return interpolate_age(0, age);
  if (Mstar_init >= init_mass[N_TRACKS-1])
    return interpolate_age(N_TRACKS-1, age);

  for (int k = 0; k < N_TRACKS-1; k++) 
  {
    double m0 = init_mass[k];
    double m1 = init_mass[k+1];
      if (Mstar_init >= m0 && Mstar_init <= m1) 
      {
        double y0 = interpolate_age(k,   age);
        double y1 = interpolate_age(k+1, age);
        return linear_interpolation(Mstar_init, m0, m1, y0, y1);
      }
  }
  // fallback
  return 0.0; 
}