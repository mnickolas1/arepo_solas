#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"


double evaluate_mu(int i)
{
#if defined(PRIMORDIAL_COOLING)
  return primordial_mu(i);

#elif defined(USE_GRACKLE)
  return grackle_mu(i);

#else
  /* Neutral primordial gas */
  return 4.0 / (1.0 + 3.0 * HYDROGEN_MASSFRAC);
#endif
}

double evaluate_gamma(int i)
{
#if defined(PRIMORDIAL_COOLING)
  return primordial_gamma(i);

#elif defined(USE_GRACKLE)
  return grackle_gamma(i);

#else
  /* Fallback */
  return GAMMA;
#endif
}

void update_mu_gamma(void)
{ 
  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      SphP[i].Mu = evaluate_mu(i);
      SphP[i].Gamma = evaluate_gamma(i);
    }
}

double evaluate_temp(int i)
{
  double temp = (SphP[i].Utherm * All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s) 
              * SphP[i].Mu * PROTONMASS * (SphP[i].Gamma - 1.0) / BOLTZMANN;
  
  return temp;
}

double evaluate_numberdens(int i)
{
  double number_dens = (SphP[i].Density * All.cf_UnitDensity_in_cgs) / SphP[i].Mu / PROTONMASS;

  return number_dens;
}

double evaluate_pressure(int i) 
{
  double pressure = (SphP[i].Gamma - 1.0) * SphP[i].Density * SphP[i].Utherm; 
  
  return pressure;
}
