/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/cooling/cooling.c
 * \date        05/2018
 * \brief       Module for gas radiative cooling
 * \details     contains functions:
 *                double DoCooling(double u_old, double rho, double dt, double
 *                  *ne_guess)
 *                double GetCoolingTime(double u_old, double rho, double
 *                  *ne_guess)
 *                double convert_u_to_temp(double u, double rho, double
 *                  *ne_guess)
 *                void find_abundances_and_rates(double logT, double rho,
 *                  double *ne_guess)
 *                double CoolingRateFromU(double u, double rho, double
 *                  *ne_guess)
 *                void SetOutputGasState(int i, double *ne_guess, double *nH0,
 *                  double *coolrate)
 *                double CoolingRate(double logT, double rho, double *nelec)
 *                void MakeRateTable(void)
 *                void ReadIonizeParams(char *fname, int which)
 *                void IonizeParamsUVB(void)
 *                void SetZeroIonization(void)
 *                void IonizeParams(void)
 *                void InitCool(void)
 *                void cooling_only(void)
 *                void cool_cell(int i)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

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

double evaluate_temp(int i)
{
  double mu = evaluate_mu(i);
  double temp = (SphP[i].Utherm * All.UnitVelocity_in_cm_per_s*All.UnitVelocity_in_cm_per_s) 
              * mu * PROTONMASS * GAMMA_MINUS1 / BOLTZMANN;
  
  return temp;
}

double evaluate_numberdens(int i)
{
  double mu = evaluate_mu(i);
  double number_dens = (SphP[i].Density * All.cf_UnitDensity_in_cgs) / mu / PROTONMASS;

  return number_dens;
}

#ifdef COOLING

void InitCool(void)
{
#if defined(PRIMORDIAL_COOLING)
  InitPrimordialCooling();

#elif defined(USE_GRACKLE)
  InitGrackle();

#else
#endif
}

/*! \brief Apply the isochoric cooling to all the active gas cells.
 *
 *  \return void
 */
void cooling_only(void) /* normal cooling routine when star formation is disabled */
{
  int idx, i;

  CPU_Step[CPU_MISC] += measure_time();

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i >= 0)
        {
          if(P[i].Mass == 0 && P[i].ID == 0)
            continue; /* skip cells that have been swallowed or eliminated */

          cool_cell(i);
        }
    }
  CPU_Step[CPU_COOLINGSFR] += measure_time();
}

/*! \brief Apply the isochoric cooling to a given gas cell.
 *
 *  This function applies the normal isochoric cooling to a single gas cell.
 *  Once the cooling has been applied according to one of the cooling models
 *  implemented, the internal energy per unit mass, the total energy and the
 *  pressure of the cell are updated.
 *
 *  \param[in] i Index of the gas cell to which cooling is applied.
 *
 *  \return void
 */
void cool_cell(int i)
{
  double dt, dtime, dtcool;

  dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

  dtime = All.cf_atime * dt / All.cf_time_hubble_a;

  dtcool = dtime;

  double unew;

#if defined(PRIMORDIAL_COOLING)
  double dens = SphP[i].Density;
  double u_old = dmax(All.MinEgySpec, SphP[i].Utherm);
  double ne = SphP[i].Ne;

  unew = (dtcool > 0) ? DoPrimordialCooling(u_old, dens * All.cf_a3inv, dtcool, &ne) : u_old;
  
  SphP[i].Ne = ne;

#elif defined(USE_GRACKLE)
  unew = (dtcool > 0) ? CallGrackle(i, dtcool, 0) : dmax(All.MinEgySpec, SphP[i].Utherm);

#else
#endif

  if(unew < 0)
    terminate("invalid temperature: Thistask=%d i=%d unew=%g\n", ThisTask, i, unew);

  double du = unew - SphP[i].Utherm;

  if(unew < All.MinEgySpec)
    du = All.MinEgySpec - SphP[i].Utherm;

  SphP[i].Utherm += du;
  SphP[i].Energy += All.cf_atime * All.cf_atime * du * P[i].Mass;

#ifdef OUTPUT_COOLHEAT
  if(dtime > 0)
    SphP[i].CoolHeat = du * P[i].Mass / dtime;
#endif /* #ifdef OUTPUT_COOLHEAT */

  set_pressure_of_cell(i);
}

double GetCoolingTime(int i)
{
  double tcool;

#if defined(PRIMORDIAL_COOLING)
  double dens = SphP[i].Density;
  double u_old = dmax(All.MinEgySpec, SphP[i].Utherm);
  double ne = SphP[i].Ne;

  tcool = GetPrimordialCoolingTime(u_old, dens * All.cf_a3inv, &ne);

#elif defined(USE_GRACKLE)
  double LambdaNet = CallGrackle(i, 0.0, 1);

  tcool = (LambdaNet >= 0) ? 0.0 : -LambdaNet;
  
#else
#endif

  return tcool;
}

#endif