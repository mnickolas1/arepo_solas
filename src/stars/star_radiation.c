#include "../main/allvars.h"
#include "../main/proto.h"

#include "../extern/chealpix.h"

/* clang-format off */
/* Effective attenuation kappa_ext*(1 - a*<g>) [cm^2/g gas, solar Z]
   Band-averaged over Draine 2003 (renorm. WD01) MW R_V=3.1 model,
   kext_albedo_WD_MW_3.1_60_D03.all, energy and photon-weighted 4e4 K BB
   Gas mass per H = 2.311e-24 g (M_dust/H = 1.398e-26, M_gas/M_dust = 165.3) */
double Kappa_E[WAVEBANDS] = {
  [INFRARED] = 34.9,
  [OPTICAL] = 278.3,
  [ULTRAVIOLET] = 417.7,
  [LYMAN_WERNER] = 736.6, 
  [IONIZING_HI] = 899.5,
  [IONIZING_H2] = 903.5,
  [IONIZING_HeI] = 460.0,
  [IONIZING_HeII] = 256.4,
};

double Kappa_N[WAVEBANDS] = {
  [INFRARED] = 30.0,
  [OPTICAL] = 242.3,
  [ULTRAVIOLET] = 406.9,
  [LYMAN_WERNER] = 731.4, 
  [IONIZING_HI] = 898.8,
  [IONIZING_H2] = 925.8,
  [IONIZING_HeI] = 469.2,
  [IONIZING_HeII] = 257.4,
};

/* Fraction of kappa_eff-attenuated energy that is truly absorbed (heats grains):
   f_abs = kappa_abs/kappa_eff = (1-a)/(1-a<g>), D03 MW dust, band-averaged
   Remainder is non-forward-scattered light, is removed from the ray and 
   delivers momentum but does not contribute to heating */
double AbsorbedFraction[WAVEBANDS] = {
  [INFRARED] = 0.54,
  [OPTICAL] = 0.62,
  [ULTRAVIOLET] = 0.81,
  [LYMAN_WERNER] = 0.88,
  [IONIZING_HI] = 0.91,
  [IONIZING_H2] = 0.93,
  [IONIZING_HeI] = 0.94,
  [IONIZING_HeII] = 0.97,
};

/* Correction to kappa_eff-attenuated energy to express the true momentum transfer 
   kappa_ext*(1 - a*<g>) is computed with <g> = max(0, <g>) to get the correct opacity
   so undestimates backward scattered momentum tranfer 
   Not needed for this dust model */
double MomentumFraction[WAVEBANDS] = {
  [INFRARED] = 1.00,
  [OPTICAL] = 1.00,
  [ULTRAVIOLET] = 1.00,
  [LYMAN_WERNER] = 1.00,
  [IONIZING_HI] = 1.00,
  [IONIZING_HeI] = 1.00,
  [IONIZING_HeII] = 1.00,
};

/* f_rerad = f_abs*(1-eps_pe); eps_pe = 0.05 for the two UV bands only */
double ReradiatedFraction[WAVEBANDS] = {
  [INFRARED] = 0.54,
  [OPTICAL] = 0.62,
  [ULTRAVIOLET] = 0.77,
  [LYMAN_WERNER] = 0.84,
  [IONIZING_HI] = 0.91,
  [IONIZING_H2] = 0.93,
  [IONIZING_HeI] = 0.94,
  [IONIZING_HeII] = 0.97,
};

double SigmaH2 = SIGMA_DISS / F_DISS;

double Sigma_E[WAVEBANDS][N_ION_SPECIES] = {
  [INFRARED] = {0.0, 0.0, 0.0, 0.0},
  [OPTICAL] = {0.0, 0.0, 0.0, 0.0},
  [ULTRAVIOLET] = {0.0, 0.0, 0.0, 0.0},
  [LYMAN_WERNER] = {0.0, 0.0, 0.0, 0.0},
  [IONIZING_HI] = {5.3851e-18, 0.0000e+00, 0.0000e+00, 0.0000e+00},
  [IONIZING_H2] = {2.7415e-18, 6.3725e-18, 0.0000e+00, 0.0000e+00},
  [IONIZING_HeI] = {8.1308e-19, 2.9751e-18, 5.7225e-18, 0.0000e+00},
  [IONIZING_HeII] = {1.0127e-19, 2.6967e-19, 1.4614e-18, 1.3294e-18},
};

double Sigma_N[WAVEBANDS][N_ION_SPECIES] = {
  [INFRARED] = {0.0, 0.0, 0.0, 0.0},
  [OPTICAL] = {0.0, 0.0, 0.0, 0.0},
  [ULTRAVIOLET] = {0.0, 0.0, 0.0, 0.0},
  [LYMAN_WERNER] = {0.0, 0.0, 0.0, 0.0},
  [IONIZING_HI] = {5.4042e-18, 0.0000e+00, 0.0000e+00, 0.0000e+00},
  [IONIZING_H2] = {2.8600e-18, 6.2862e-18, 0.0000e+00, 0.0000e+00},
  [IONIZING_HeI] = {8.4911e-19, 3.1108e-18, 5.8894e-18, 0.0000e+00},
  [IONIZING_HeII] = {1.0236e-19, 2.7276e-19, 1.4735e-18, 1.3425e-18},
};
/* clang-format on */

/* WG19 self-shielding exponent, inside the fit range */
static inline double h2shield_alpha(double temp, double n)
{
  const double lT = log10(fmin(fmax(temp, H2_SHIELD_TMIN), H2_SHIELD_TMAX));
  const double ln = log10(fmin(fmax(n, H2_SHIELD_NMIN), H2_SHIELD_NMAX));

  const double alpha = (0.8711 * lT - 1.928) * exp(-0.2856 * ln) + (-0.9639 * lT + 3.892);

  return fmax(alpha, 0.0);
}

/* H2 thermal Doppler parameter in km/s, b = sqrt(2 k T / m_H2); T floored to keep b > 0 */
static inline double h2shield_b5(double temp)
{
  return 1.0e-5 * sqrt(BOLTZMANN * fmax(temp, 1.0) / PROTONMASS);
}

void update_opac(void)
{
  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      double Units;

      Units = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm / All.cf_UnitMass_in_g;

      double Density = (P[i].Mass + SphP[i].StarMassFeed) / SphP[i].Volume;

#ifdef METALS
      double Zsol = ((SphP[i].GasMetals + SphP[i].StarMetalsFeed) / (P[i].Mass + SphP[i].StarMassFeed)) / SOLAR_METALLICITY;
      SphP[i].OpacityScaling[CH_DUST] = dust_to_gas_ratio(Zsol) / DUST_TO_GAS_RATIO * Density / Units;
#else
      SphP[i].OpacityScaling[CH_DUST] = 0;
#endif

      Units = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm;

      double n_H2 = SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / SphP[i].Volume / (2 * PROTONMASS / All.cf_UnitMass_in_g);

      SphP[i].OpacityScaling[CH_LWH2] = fmax(0.0, n_H2 / Units);

      /* Shielding parameters for the local gas */
      double temp = evaluate_temp(i);
      double number_dens = evaluate_numberdens(i);

      SphP[i].H2ShieldAlpha = h2shield_alpha(temp, number_dens);
      SphP[i].H2ShieldB5 = h2shield_b5(temp);

      double n_Ionizing[4] = {
          SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / SphP[i].Volume / (PROTONMASS / All.cf_UnitMass_in_g),
          SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / SphP[i].Volume / (2 * PROTONMASS / All.cf_UnitMass_in_g),
          SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / SphP[i].Volume / (4 * PROTONMASS / All.cf_UnitMass_in_g),
          SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / SphP[i].Volume / (4 * PROTONMASS / All.cf_UnitMass_in_g)};

      for(int s = 0; s < N_ION_SPECIES; s++)
        SphP[i].OpacityScaling[CH_HI + s] = fmax(0.0, n_Ionizing[s] / Units);
    }
}

#ifdef IR_MOMENTUM_BOOST
double dtau_IR(int i, double length)
{
  double kappa_rerad = 1.0;

  double Dtau_IR = All.IRDtauMomentumBoostCoeff * kappa_rerad * SphP[i].OpacityScaling[CH_DUST] * length;

  return Dtau_IR;
}
#endif

/* (exp(s*y) - 1) / s, continuous through s = 0 */
static inline double expm1_over(double s, double y)
{
  const double sy = s * y;
  return fabs(sy) == 0 ? y : expm1(sy) / s;
}

/* Fraction of the LW band absorbed in H2 lines between N_H2 and N_H2 + dN_H2,
   dA = SigmaH2 * int f_sh(N'; alpha, b5) dN', in closed form
   Written in differences so dN << N does not cancel */
double h2shield_dA(double N_H2, double dN_H2, double alpha, double b5)
{
  if(dN_H2 <= 0.0)
    return 0.0;

  const double x1 = N_H2 / H2_SHIELD_N0;
  const double dx = dN_H2 / H2_SHIELD_N0;

  /* 0.965 / (1 + x/b5)^alpha */
  const double s = 1.0 - alpha;
  const double t1 = 0.965 * b5 * exp(s * log1p(x1 / b5)) * expm1_over(s, log1p(dx / (b5 + x1)));

  /* 0.035 / sqrt(1 + x) * exp(-a sqrt(1 + x)) */
  const double u1 = sqrt(1.0 + x1);
  const double u2 = sqrt(1.0 + x1 + dx);
  const double t2 = -(0.07 / H2_SHIELD_A) * exp(-H2_SHIELD_A * u1) * expm1(-H2_SHIELD_A * dx / (u1 + u2));

  return SigmaH2 * H2_SHIELD_N0 * (t1 + t2);
}

/* Effective optical depth of a cell for H2 lines: exp(-dtau) = (1 - A_H2 - dA) / (1 - A_H2) */
double h2shield_dtau(double A_H2, double dA)
{
  const double T = 1.0 - A_H2;

  if(dA <= 0.0)
    return 0.0;

  if(dA >= T)
    return RAD_TAU_SAT;

  return fmin(-log1p(-dA / T), RAD_TAU_SAT);
}

/* Helpers for rotation */
static inline unsigned long long splitmix64(unsigned long long *s)
{
  unsigned long long z = (*s += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

static inline double splitmix_uniform(unsigned long long *s)
{
  return (splitmix64(s) >> 11) * (1.0 / 9007199254740992.0); /* [0,1) */
}

/* Shoemake (1992): uniform random quaternion -> rotation matrix */
static void get_ray_rotation(unsigned long long seed, double R[3][3])
{
  unsigned long long s = seed * 0x2545F4914F6CDD1Dull + 0x9E3779B97F4A7C15ull;

  double u1 = splitmix_uniform(&s);
  double u2 = splitmix_uniform(&s);
  double u3 = splitmix_uniform(&s);

  double r1 = sqrt(1.0 - u1), r2 = sqrt(u1);
  double t1 = 2.0 * M_PI * u2, t2 = 2.0 * M_PI * u3;

  double x = r1 * sin(t1), y = r1 * cos(t1);
  double z = r2 * sin(t2), w = r2 * cos(t2);

  R[0][0] = 1.0 - 2.0 * (y * y + z * z);
  R[0][1] = 2.0 * (x * y - z * w);
  R[0][2] = 2.0 * (x * z + y * w);

  R[1][0] = 2.0 * (x * y + z * w);
  R[1][1] = 1.0 - 2.0 * (x * x + z * z);
  R[1][2] = 2.0 * (y * z - x * w);

  R[2][0] = 2.0 * (x * z - y * w);
  R[2][1] = 2.0 * (y * z + x * w);
  R[2][2] = 1.0 - 2.0 * (x * x + y * y);
}

/* Seed from the star's global ID: stable under domain decomposition and
   identical on every rank, so a ray splitting on an imported domain
   reproduces the same children as its home rank would have */
static inline unsigned long long seed_rotation(MyIDType id)
{
  unsigned long long seed = (unsigned long long)id;
  seed = seed * 0x9E3779B97F4A7C15ull + (unsigned long long)All.Ti_Current;
  return seed;
}

/* Rotated HEALPix direction */
static inline void healpix_dir(unsigned long long seed, int nside, int ipix, double *dir)
{
  static unsigned long long cached_seed = 0;
  static int cached_valid = 0;
  static double R[3][3];

  if(!cached_valid || seed != cached_seed)
    {
      get_ray_rotation(seed, R);
      cached_seed = seed;
      cached_valid = 1;
    }

  double v[3];
  pix2vec_nest(nside, ipix, v);

  for(int k = 0; k < 3; k++)
    dir[k] = R[k][0] * v[0] + R[k][1] * v[1] + R[k][2] * v[2];
}

static RayWorkStack *init_work_stack(long long capacity)
{
  RayWorkStack *w = malloc(sizeof(RayWorkStack));

  w->n = 0;
  w->capacity = capacity;
  w->rays = malloc(capacity * sizeof(RayPacket));
  return w;
}

void append_ray(RayWorkStack *w, const RayPacket *ray)
{
  if(w->n >= w->capacity)
    {
      w->capacity *= 2;
      w->rays = realloc(w->rays, w->capacity * sizeof(RayPacket));
    }
  w->rays[w->n++] = *ray;
}

static void free_work_stack(RayWorkStack *w)
{
  free(w->rays);
  free(w);
}

/*
 * Rays are born inside their host cell, which the star density
 * machinery already identified through the SPH neighbour loop
 */
static void init_rays(RayWorkStack *work)
{
  double SQRT3 = sqrt(3.0);

  for(int ev = 0; ev < MechanicalFeedbackEvents.NumEvents;)
    {
      int host = MechanicalFeedbackEvents.MechanicalFeedbackData[ev].HostIndex;

#ifdef STAR_IN_CELL

      /* Superpose every star in this cell into one source on the host generator */
      WavebandData Radiated_Cell[WAVEBANDS];

      for(int w = 0; w < WAVEBANDS; w++)
        Radiated_Cell[w].Energy = Radiated_Cell[w].Photons = 0.0;

      for(int h = 0; h < SphP[host].Host; h++)
        {
          Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h].MechanicalFeedback;

          for(int w = 0; w < WAVEBANDS; w++)
            {
              Radiated_Cell[w].Energy += MechanicalFeedback->Radiated[w].Energy;
              Radiated_Cell[w].Photons += MechanicalFeedback->Radiated[w].Photons;
            }
        }

      /* Skip dark stars entirely rather than pushing dead rays */
      int flag_luminosity = 0;
      for(int w = 0; w < WAVEBANDS; w++)
        {
          if(Radiated_Cell[w].Energy > 0.0 || Radiated_Cell[w].Photons > 0.0)
            {
              flag_luminosity = 1;
              break;
            }
        }

      if(flag_luminosity)
        {
          /* Loop over rays for this host */
          for(int iray = 0; iray < NRays; iray++)
            {
              RayPacket ray = {0};

              ray.star_id = P[host].ID;

              ray.cell = host;

              ray.pos[0] = 0.0;
              ray.pos[1] = 0.0;
              ray.pos[2] = 0.0;

              unsigned long long rotation_seed = seed_rotation(ray.star_id);
              healpix_dir(rotation_seed, NSIDE_MIN, iray, ray.dir);

              ray.t = 0.0;
              ray.t_maximum = All.RayMaxDistance > 0 ? All.RayMaxDistance : SQRT3 * All.BoxSize;

              ray.nside = NSIDE_MIN;
              ray.healpix_pixel = iray;

              ray.active_bands = NO_IR_ACTIVE;

              for(int w = 0; w < WAVEBANDS; w++)
                {
                  ray.Radiated[w].Energy = Radiated_Cell[w].Energy / NRays;
                  ray.Radiated[w].Photons = Radiated_Cell[w].Photons / NRays;

#ifndef RAD_TOTAL_TRUNCATION
                  ray.Radiated_Init[w].Energy = Radiated_Cell[w].Energy / NRays;
                  ray.Radiated_Init[w].Photons = Radiated_Cell[w].Photons / NRays;
#endif

                  if(ray.Radiated[w].Energy <= 0.0 && ray.Radiated[w].Photons <= 0.0)
                    ray.active_bands &= (uint8_t)(~(1u << w));
                }

#ifdef RAD_TOTAL_TRUNCATION
              ray.E_init = ray.N_init = 0.0;
              for(int w = 0; w < WAVEBANDS; w++)
                {
                  if(!(ray.active_bands & (1u << w)))
                    continue;

                  ray.E_init += ray.Radiated[w].Energy;

                  if((BandTrackPhotons >> w) & 1u)
                    ray.N_init += ray.Radiated[w].Photons;
                }
#endif

              ray.N_H2 = 0.0;
              ray.A_H2 = 0.0;

              if(ray.active_bands == 0)
                continue;

#ifdef RT_STATISTICS
              rt_statistics_init(&ray);
#endif

              append_ray(work, &ray);
            }
        }

#else

      double xtmp, ytmp, ztmp;

      for(int h = 0; h < SphP[host].Host; h++)
        {
          Mechanical_Feedback_Data *MechanicalFeedbackData = &MechanicalFeedbackEvents.MechanicalFeedbackData[ev + h];
          Mechanical_Feedback *MechanicalFeedback = &MechanicalFeedbackData->MechanicalFeedback;

          /* Skip dark stars entirely rather than pushing dead rays */
          int flag_luminosity = 0;
          for(int w = 0; w < WAVEBANDS; w++)
            {
              if(MechanicalFeedback->Radiated[w].Energy > 0.0 || MechanicalFeedback->Radiated[w].Photons > 0.0)
                {
                  flag_luminosity = 1;
                  break;
                }
            }

          if(!flag_luminosity)
            continue;

          /* Star position relative to the host generator, minimum image */
          double xrel[3];

          xrel[0] = NEAREST_X(MechanicalFeedback->StarPosition[0] - P[host].Pos[0]);
          xrel[1] = NEAREST_Y(MechanicalFeedback->StarPosition[1] - P[host].Pos[1]);
          xrel[2] = NEAREST_Z(MechanicalFeedback->StarPosition[2] - P[host].Pos[2]);

          /* Loop over rays for this star */
          for(int iray = 0; iray < NRays; iray++)
            {
              RayPacket ray = {0};

              ray.star_id = MechanicalFeedbackData->StarParticleID;

              ray.cell = host;

              ray.pos[0] = xrel[0];
              ray.pos[1] = xrel[1];
              ray.pos[2] = xrel[2];

              unsigned long long rotation_seed = seed_rotation(ray.star_id);
              healpix_dir(rotation_seed, NSIDE_MIN, iray, ray.dir);

              ray.t = 0.0;
              ray.t_maximum = All.RayMaxDistance > 0 ? All.RayMaxDistance : SQRT3 * All.BoxSize;

              ray.nside = NSIDE_MIN;
              ray.healpix_pixel = iray;

              ray.active_bands = NO_IR_ACTIVE;

              for(int w = 0; w < WAVEBANDS; w++)
                {
                  ray.Radiated[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
                  ray.Radiated[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;

#ifndef RAD_TOTAL_TRUNCATION
                  ray.Radiated_Init[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
                  ray.Radiated_Init[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;
#endif

                  if(ray.Radiated[w].Energy <= 0.0 && ray.Radiated[w].Photons <= 0.0)
                    ray.active_bands &= (uint8_t)(~(1u << w));
                }

#ifdef RAD_TOTAL_TRUNCATION
              ray.E_init = ray.N_init = 0.0;
              for(int w = 0; w < WAVEBANDS; w++)
                {
                  if(!(ray.active_bands & (1u << w)))
                    continue;

                  ray.E_init += ray.Radiated[w].Energy;

                  if((BandTrackPhotons >> w) & 1u)
                    ray.N_init += ray.Radiated[w].Photons;
                }
#endif

              ray.N_H2 = 0.0;
              ray.A_H2 = 0.0;

              if(ray.active_bands == 0)
                continue;

#ifdef RT_STATISTICS
              rt_statistics_init(&ray);
#endif

              append_ray(work, &ray);
            }
        }

#endif

      ev += SphP[host].Host;
    }
}

/* Splits to 4 child rays
   Children inherit position, cell and path length,
   so they simply restart the exit search in the cell the parent entered */
void split_ray(const RayPacket *parent, RayPacket children[4])
{
  int new_nside = parent->nside * 2;

  for(int k = 0; k < 4; k++)
    {
      /* Copy all state including t, cell, pos, N_H2, active_bands */
      children[k] = *parent;

      children[k].nside = new_nside;
      children[k].healpix_pixel = 4 * parent->healpix_pixel + k;
      children[k].locate_head = 1;

      unsigned long long rotation_seed = seed_rotation(parent->star_id);
      healpix_dir(rotation_seed, new_nside, children[k].healpix_pixel, children[k].dir);

      for(int w = 0; w < WAVEBANDS; w++)
        {
          children[k].Radiated[w].Energy = parent->Radiated[w].Energy * 0.25;
          children[k].Radiated[w].Photons = parent->Radiated[w].Photons * 0.25;

#ifndef RAD_TOTAL_TRUNCATION
          children[k].Radiated_Init[w].Energy = parent->Radiated_Init[w].Energy * 0.25;
          children[k].Radiated_Init[w].Photons = parent->Radiated_Init[w].Photons * 0.25;
#endif
        }

#ifdef RAD_TOTAL_TRUNCATION
      children[k].E_init = parent->E_init * 0.25;
      children[k].N_init = parent->N_init * 0.25;
#endif
    }

#ifdef RT_STATISTICS
  RTStatisticsLocal.n_split++;
#endif
}

/* Sparse, neighbour-restricted ray exchange */
int RayNgbNTask = 0;
int *RayNgbTask = NULL;
int *RayTaskToNgb = NULL;

/*
 * Mesh-neighbour graph - shared by both back ends
 *
 * Walk every local cell's Delaunay connection list
 * and flag the ranks that own a face-defining neighbour
 */
void ray_neighbours_init(void)
{
  char *sflag = malloc(NTask * sizeof(char));
  char *rflag = malloc(NTask * sizeof(char));

  memset(sflag, 0, NTask * sizeof(char));

  for(int i = 0; i < NumGas; i++)
    {
      int q = SphP[i].first_connection;

      while(q >= 0)
        {
          if(q >= MaxNvc)
            terminate("ray_neighbours_init(): strange connectivity q=%d MaxNvc=%d cell=%d!\n", q, MaxNvc, i);

          const int dp = DC[q].dp_index;

          if(Mesh.DP[dp].index >= 0)
            {
              const int t = DC[q].task;

              if(t < 0 || t >= NTask)
                terminate("ray_neighbours_init(): DC[%d].task = %d out of range on cell%d!\n", q, t, i);

              if(t != ThisTask)
                sflag[t] = 1;
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }
    }

  /*
   * Symmetrise: this rank must be able to RECEIVE from anyone who can send to it
   * (AREPO's face connectivity should already be symmetric so sflag == rflag)
   */
  MPI_Alltoall(sflag, 1, MPI_CHAR, rflag, 1, MPI_CHAR, MPI_COMM_WORLD);

  RayNgbNTask = 0;
  for(int t = 0; t < NTask; t++)
    {
      if(sflag[t] || rflag[t])
        RayNgbNTask++;
    }

  RayNgbTask = malloc((RayNgbNTask > 0 ? RayNgbNTask : 1) * sizeof(int));
  RayTaskToNgb = malloc(NTask * sizeof(int));

  for(int t = 0; t < NTask; t++)
    RayTaskToNgb[t] = -1;

  int k = 0;
  for(int t = 0; t < NTask; t++)
    {
      if(sflag[t] || rflag[t])
        {
          RayTaskToNgb[t] = k;
          RayNgbTask[k++] = t;
        }
    }

  free(rflag);
  free(sflag);

  int ngb_max, ngb_sum;
  MPI_Allreduce(&RayNgbNTask, &ngb_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&RayNgbNTask, &ngb_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  mpi_printf("STAR_RADIATION: RayPacket = %d B, comm neighbours: mean %d, max %d (of %d ranks)\n", (int)sizeof(RayPacket),
             ngb_sum / NTask, ngb_max, NTask);
}

void ray_neighbours_free(void)
{
  free(RayTaskToNgb);
  RayTaskToNgb = NULL;
  free(RayNgbTask);
  RayNgbTask = NULL;
  RayNgbNTask = 0;
}

static void radiation_feedback(void)
{
  /* Indexed by ionizing species a = 0,1,2,3 (HI, H2, HeI, HeII) */
  static const double IonThreshold_eV[N_ION_SPECIES] = {13.6, 15.4, 24.6, 54.4};
  static const int IonGrackle[N_ION_SPECIES] = {GRACKLE_HI, GRACKLE_H2I, GRACKLE_HeI, GRACKLE_HeII};
  static const double IonAtomicMass[N_ION_SPECIES] = {1.0, 2.0, 4.0, 4.0};

  const double L3 = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm;

  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      const double V = SphP[i].Volume;
      const double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      if(dt <= 0.0 || V <= 0.0)
        goto reset;

      const double V_cgs = V * L3;
      const double dt_cgs = dt * All.cf_UnitTime_in_s;

#ifdef PHOTOELECTRIC_HEATING
      const double epsilon_pe = 0.05;
      const double E_pe = SphP[i].AbsorbedPE * epsilon_pe * All.cf_UnitEnergy_in_cgs;

      SphP[i].PE_VolHeatingRate += E_pe / dt_cgs / V_cgs;
#endif

#ifdef DISSOCIATION
      const double n_H2 = SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / V / (2.0 * PROTONMASS / All.cf_UnitMass_in_g);

      /* AbsorbedH2Line holds pumped photons; F_DISS is the branching ratio */
      if(n_H2 > 0.0)
        SphP[i].H2_DissociationRate += F_DISS * SphP[i].AbsorbedH2Line / (dt / All.cf_hubble_a) / V / n_H2;
#endif

#ifdef PHOTOIONIZATION
      for(int s = 0; s < N_ION_SPECIES; s++)
        {
          const double n = SphP[i].GrackleSpeciesConserved(IonGrackle[s]) / V / (IonAtomicMass[s] * PROTONMASS / All.cf_UnitMass_in_g);

          const double N_abs = SphP[i].AbsorbedIonizing[s].Photons;
          const double E_abs = SphP[i].AbsorbedIonizing[s].Energy * All.cf_UnitEnergy_in_cgs;

          const double E_exc = E_abs - N_abs * IonThreshold_eV[s] * ELECTRONVOLT_IN_ERGS;

          if(n <= 0.0)
            continue;

          const double n_cgs = n / L3;

          if(E_exc > 0.0)
            SphP[i].IonHeatingRate[s] += E_exc / dt_cgs / V_cgs / n_cgs;
          else if(N_abs > 0.0)
            warn(
                "STAR_RADIATION: sub-threshold mean photon energy, species %d, cell %d "
                "(E_abs=%g N_abs=%g) \n",
                s, i, E_abs, N_abs);

          SphP[i].IonizationRate[s] += N_abs / (dt / All.cf_hubble_a) / V / n;
        }
#endif

    reset:

#ifdef PHOTOELECTRIC_HEATING
      SphP[i].AbsorbedPE = 0.0;
#endif

#ifdef DISSOCIATION
      SphP[i].AbsorbedH2Line = 0.0;
#endif

#ifdef PHOTOIONIZATION
      for(int s = 0; s < N_ION_SPECIES; s++)
        SphP[i].AbsorbedIonizing[s].Energy = SphP[i].AbsorbedIonizing[s].Photons = 0.0;
#endif
    }
}

#ifdef RT_TIMESTEP
static void rt_timestep(void)
{
  const double eps_ion = All.RTIonizationTimestepFraction;

  for(int idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      int i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* Total hydrogen */
      double m_H = SphP[i].GrackleSpeciesConserved(GRACKLE_HI)
                 + SphP[i].GrackleSpeciesConserved(GRACKLE_HII)
                 + SphP[i].GrackleSpeciesConserved(GRACKLE_H2I)
                 + SphP[i].GrackleSpeciesConserved(GRACKLE_H2II)
                 + SphP[i].GrackleSpeciesConserved(GRACKLE_HM);

      /* Total helium */
      double m_He = SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) 
                  + SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) 
                  + SphP[i].GrackleSpeciesConserved(GRACKLE_HeIII);

      double rate = 0.0;

      if(m_H > 0)
        {
#ifdef PHOTOIONIZATION
          double x_HI = SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / m_H;
          rate = fmax(rate, SphP[i].IonizationRate[SP_HI] * x_HI);
#endif

          double rate_H2 = 0.0;

#ifdef DISSOCIATION
          rate_H2 += SphP[i].H2_DissociationRate;
#endif

#ifdef PHOTOIONIZATION
          rate_H2 += SphP[i].IonizationRate[SP_H2];
#endif
          double x_H2 = SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / m_H;
          rate = fmax(rate, rate_H2 * x_H2);
        }

      if(m_He > 0)
        {
#ifdef PHOTOIONIZATION
          double x_HeI = SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / m_He;
          double x_HeII = SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / m_He;
 
          rate = fmax(rate, SphP[i].IonizationRate[SP_HeI] * x_HeI);
          rate = fmax(rate, SphP[i].IonizationRate[SP_HeII] * x_HeII);
#endif
        }


      SphP[i].RT_Timestep = (rate > 0.0) ? eps_ion / rate : All.MaxSizeTimestep / All.cf_hubble_a;
    }
}
#endif

void star_radiation(void)
{
  TIMER_START(CPU_STARS_RADIATION);

  double t0, t1;

  update_opac();

  /* Zero accumulators before the walk */
  for(int i = 0; i < NumGas; i++)
    {
#ifdef PHOTOELECTRIC_HEATING
      SphP[i].AbsorbedPE = 0.0;
#endif

#ifdef DISSOCIATION
      SphP[i].AbsorbedH2Line = 0.0;
#endif

#ifdef PHOTOIONIZATION
      for(int s = 0; s < N_ION_SPECIES; s++)
        SphP[i].AbsorbedIonizing[s].Energy = SphP[i].AbsorbedIonizing[s].Photons = 0.0;
#endif
    }

  long long n_sources_local = 0;

#ifdef STAR_IN_CELL
  for(int ev = 0; ev < MechanicalFeedbackEvents.NumEvents;)
    {
      int host = MechanicalFeedbackEvents.MechanicalFeedbackData[ev].HostIndex;
      n_sources_local++;
      ev += SphP[host].Host;
    }
#else
  n_sources_local = MechanicalFeedbackEvents.NumEvents;
#endif

  long long n_rays_local = n_sources_local * NRays;

  long long n_rays_global;
  sumup_longs(1, &n_rays_local, &n_rays_global);

  mpi_printf("STAR_RADIATION: Initializing radiation with %12lld rays\n", n_rays_global);

  /* Floor so ranks with no local stars still have a buffer to receive imports */
  long long work_capacity = n_rays_local > 0 ? 4 * n_rays_local : 1024;

  RayWorkStack *work = init_work_stack(work_capacity);
  RayComms *comm = ray_comms_init(work);

#ifdef RT_STATISTICS
  rt_statistics_reset();
#endif

  init_rays(work);

  t0 = second();

  ray_comms_walk(work, comm);

  t1 = second();
  mpi_printf("STAR_RADIATION: walk complete (%g sec)\n", timediff(t0, t1));

#ifdef RT_STATISTICS
  rt_statistics_report();
#endif

  ray_comms_free(comm);
  free_work_stack(work);

  radiation_feedback();

#ifdef RT_TIMESTEP
  rt_timestep();
#endif

  TIMER_STOP(CPU_STARS_RADIATION);
}