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
 *              solar metallicity.  Shared by the radiative transfer dust
 *              opacities (stars/star_radiation.c) and by the dust density
 *              handed to grackle (cooling/grackle.c), so that the two cannot
 *              disagree about how much dust a cell holds.
 *
 *              contains functions:
 *                void init_dust_to_gas(void)
 *                double dust_to_gas_factor(double Zsol)
 *                double Zsol_from_dust_to_gas_factor(double factor)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */

#include <math.h>
#include <stdio.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#ifdef METALS

/*! \brief Metallicity in solar units at the break of the broken power law.
 *
 *  Z_t/Z_sun = 10^(x_t - x_sun) with x = 12 + log(O/H).  Cached by
 *  init_dust_to_gas() so that the pow() is not repeated for every cell.
 */
static double ZsolBreak;

/*! \brief The dust-to-gas factor at the break, ZsolBreak^DustToGasSlopeHigh. */
static double FactorBreak;

/*! \brief Precompute the break of the dust-to-gas relation and report it.
 *
 *  Must be called after the parameter file has been read and before the first
 *  call to dust_to_gas_factor().
 *
 *  \return void
 */
void init_dust_to_gas(void)
{
  ZsolBreak = pow(10.0, All.DustToGasBreakOH - All.DustToGasSolarOH);
  FactorBreak = pow(ZsolBreak, All.DustToGasSlopeHigh);

  /* Report the relation, and the round trip through the inverse, so that a mis-read parameter file is
   * obvious in the log.  The round trip is exact except where MinDustToGasFactor clips the result */

  const double Ztable[] = {1.0, 0.5, ZsolBreak, 0.1, 0.01};

  for(int k = 0; k < (int)(sizeof(Ztable) / sizeof(Ztable[0])); k++)
    {
      double factor = dust_to_gas_factor(Ztable[k]);
    }
}

/*! \brief Dust-to-gas mass ratio relative to its value at solar metallicity.
 *
 *  Remy-Ruyer et al. 2014 (A&A 563, A31) fit the gas-to-dust mass ratio of local galaxies as
 *  log10(G/D) = a + alpha * (x_sun - x) with x = 12 + log(O/H), and find that a single slope cannot
 *  describe the full 2 dex range: below x_t the dust content falls far faster than linearly in
 *  metallicity.  AREPO carries only a total metal mass fraction, so O/H is taken proportional to Z and
 *  x_sun - x = -log10(Z/Z_sun).  Normalising out the solar anchor a, which then cancels entirely,
 *
 *      f(Z) = (Z/Z_sun)^alpha_high                              for Z >= Z_t
 *           = (Z_t/Z_sun)^alpha_high * (Z/Z_t)^alpha_low        for Z <  Z_t
 *
 *  so f(Z_sun) = 1 by construction, and alpha_high = alpha_low = 1 recovers the constant
 *  dust-to-metals ratio that the code assumed previously.
 *
 *  Because the relation is normalised to solar, each caller divides by whatever dust-to-gas ratio its
 *  own rates or opacities were calibrated at: DUST_TO_GAS_RATIO for the Draine Kappa_E/Kappa_N tables,
 *  grackle's local_dust_to_gas_ratio for its H2-on-dust and gas-grain rates.
 *
 *  \param[in] Zsol Metallicity in solar units, Z/Z_sun.  Values <= 0, which advection undershoot does
 *             produce, return the floor rather than reaching pow() with a negative base and a
 *             fractional exponent, which would be NaN.
 *
 *  \return Dust-to-gas mass ratio divided by its value at solar metallicity.
 */
double dust_to_gas_factor(double Zsol)
{
  if(Zsol <= 0.0)
    return All.MinDustToGasFactor;

  double factor;

  if(Zsol >= ZsolBreak)
    factor = pow(Zsol, All.DustToGasSlopeHigh);
  else
    factor = FactorBreak * pow(Zsol / ZsolBreak, All.DustToGasSlopeLow);

  return fmax(All.MinDustToGasFactor, factor);
}

/*! \brief Metallicity implied by a dust-to-gas ratio; the inverse of dust_to_gas_factor().
 *
 *  The broken power law is monotonic in Z, so this inverts it exactly.  It cannot undo the
 *  MinDustToGasFactor floor: every metallicity below the floor maps to the same factor, so a factor at
 *  the floor has no unique inverse.
 *
 *  Nothing in the code calls this yet.  It is the dust -> metallicity direction, used by the
 *  verification round trip, and needed if a dust abundance is ever advected in its own right or if an
 *  effective metallicity is wanted from an observed dust-to-gas ratio.
 *
 *  \param[in] factor Dust-to-gas mass ratio relative to its value at solar metallicity.
 *
 *  \return Metallicity in solar units, Z/Z_sun.
 */
double Zsol_from_dust_to_gas_factor(double factor)
{
  if(factor <= 0.0)
    return 0.0;

  if(factor >= FactorBreak)
    return pow(factor, 1.0 / All.DustToGasSlopeHigh);

  return ZsolBreak * pow(factor / FactorBreak, 1.0 / All.DustToGasSlopeLow);
}

#endif /* #ifdef METALS */