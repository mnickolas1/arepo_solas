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
 * \file        src/cooling/cooling_proto.h
 * \date        05/2018
 * \brief       Header for cooling functions.
 * \details     Declares the model-independent gas-state helpers in cooling.c
 *              and, guarded by the corresponding config flag, the entry points
 *              of the two cooling models (primordial_cooling.c, grackle.c).
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef INLINE_FUNC
#define INLINE_FUNC
#endif /* #ifndef INLINE_FUNC */

/* --- common gas state; cooling/cooling.c is compiled unconditionally --- */

double evaluate_mu(int i);
double evaluate_temp(int i);
double evaluate_numberdens(int i);

#ifdef COOLING
void InitCool(void);

void cooling_only(void);
void cool_cell(int i);

double GetCoolingTime(int i);
#endif /* #ifdef COOLING */

/* --- primordial (Katz/KWH) network; cooling/primordial_cooling.c --- */

#ifdef PRIMORDIAL_COOLING
double primordial_mu(int i);

void InitPrimordialCooling(void);
double DoPrimordialCooling(double u_old, double rho, double dt, double *ne_guess);
double GetPrimordialCoolingTime(double u_old, double rho, double *ne_guess);

void SetOutputGasState(int i, double *ne_guess, double *nH0, double *coolrate);
void SetZeroIonization(void);
void IonizeParams(void);
#endif /* #ifdef PRIMORDIAL_COOLING */

/* --- Grackle; cooling/grackle.c --- */

#ifdef USE_GRACKLE
#include <grackle.h>

double grackle_mu(int i);

void InitGrackle(void);
double CallGrackle(int i, double u, double rho, double dt, int mode);
#endif /* #ifdef USE_GRACKLE */