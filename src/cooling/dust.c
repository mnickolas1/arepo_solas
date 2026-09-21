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
 * \file        src/cooling/dust.c
 * \date        09/2026
 * \brief       Conversion between gas metallicity and dust content.
 * \details     Implements the broken power law of Remy-Ruyer et al. 2014
 *              (A&A 563, A31) relating the gas-to-dust mass ratio to
 *              metallicity, expressed as a factor normalised to unity at
 *              solar metallicity.  
 *              Shared by the radiative transfer dust
 *              opacities (stars/star_radiation.c) and by the dust density
 *              handed to grackle (cooling/grackle.c), so that the two cannot
 *              disagree about how much dust a cell holds.
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <stdio.h>

#include "../main/allvars.h"
#include "../main/proto.h"


#define ZSOL_BREAK 0.2570395783

/* Gas to dust ratio; broken power law from Remy-Ruyer (2014) */
double dust_to_gas_ratio(double Zsol)
{
  if(Zsol <= 0.0)
    return 0.0;

  double GtoD, DtoG;

  if(Zsol >= ZSOL_BREAK)
    GtoD = 162.0 * pow(Zsol, -1.0);
  else
    GtoD = 9.12 * pow(Zsol, -3.1);

  DtoG = 1.0 / GtoD;

  return fmax(0.0, DtoG);
}