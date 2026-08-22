#include "../main/allvars.h"
#include "../main/proto.h"

#include "../extern/chealpix.h"


/* Active bands - change at init_rays */
#define ALL_BANDS_ACTIVE ((uint8_t)((1u << WAVEBANDS) - 1u))
#define NO_IR_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ~(1u << INFRARED)))
#define NO_IONIZING_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ~(1u << IONIZING_HI) & ~(1u << IONIZING_HeI) & ~(1u << IONIZING_HeII)))
#define ONLY_IONIZING_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ((1u << IONIZING_HI) | (1u << IONIZING_HeI) | (1u << IONIZING_HeII))))

/* Effective attenuation kappa_ext*(1 - a*<g>) [cm^2/g gas, solar Z]
   Band-averaged over Draine 2003 (renorm. WD01) MW R_V=3.1 model,
   kext_albedo_WD_MW_3.1_60_D03.all, energy and photon-weighted 4e4 K BB.
   Gas mass per H = 2.311e-24 g (M_dust/H = 1.398e-26, M_gas/M_dust = 165.3) */
double Kappa_E[WAVEBANDS] =
{
  [INFRARED] = 34.9,
  [OPTICAL] = 278.3,
  [ULTRAVIOLET] = 417.7,
  [LYMAN_WERNER] = 736.6, 
  [IONIZING_HI] = 902.4,
  [IONIZING_HeI] = 460.0,
  [IONIZING_HeII] = 256.4,
};

double Kappa_N[WAVEBANDS] =
{
  [INFRARED] = 30.0,
  [OPTICAL] = 242.3,
  [ULTRAVIOLET] = 406.9,
  [LYMAN_WERNER] = 731.4, 
  [IONIZING_HI] = 917.2,
  [IONIZING_HeI] = 469.2,
  [IONIZING_HeII] = 257.4,
};

/* Fraction of kappa_eff-attenuated energy that is truly absorbed (heats grains):
   f_abs = kappa_abs/kappa_eff = (1-a)/(1-a<g>), D03 MW dust, band-averaged.
   Remainder is non-forward-scattered light: removed from the ray and it does
   deliver momentum (kappa_eff is exactly the momentum-transfer opacity), but
   it must NOT contribute to heating. */
double TrueAbsorbedFraction[WAVEBANDS] =
{
  [INFRARED] = 0.54,
  [OPTICAL] = 0.62,
  [ULTRAVIOLET] = 0.81,
  [LYMAN_WERNER] = 0.88,
  [IONIZING_HI] = 0.92,
  [IONIZING_HeI] = 0.94,
  [IONIZING_HeII] = 0.97,
};

/* f_rerad = f_abs*(1-eps_pe); eps_pe = 0.05 for the two FUV bands only */
double ReradiatedFraction[WAVEBANDS] =
{
  [INFRARED] = 0.54,
  [OPTICAL] = 0.62,
  [ULTRAVIOLET] = 0.77,
  [LYMAN_WERNER] = 0.84,
  [IONIZING_HI] = 0.92,
  [IONIZING_HeI] = 0.94,
  [IONIZING_HeII] = 0.97,
};

double SigmaH2 = SIGMA_DISS / F_DISS;

double Sigma_E[3][3] = {
  { 3.4457e-18, 0.0000e+00, 0.0000e+00 },
  { 8.1308e-19, 5.7225e-18, 0.0000e+00 },
  { 1.0127e-19, 1.4614e-18, 1.3294e-18 },
};

double Sigma_N[3][3] = {
  { 3.6742e-18, 0.0000e+00, 0.0000e+00 },
  { 8.4911e-19, 5.8894e-18, 0.0000e+00 },
  { 1.0236e-19, 1.4735e-18, 1.3425e-18 },
};

void update_opac(void)
{
  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;
      
      double Units;

      Units = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm / All.cf_UnitMass_in_g;

#ifdef METALS
      double Zsol = ((SphP[i].GasMetals + SphP[i].StarMetalsFeed) / (P[i].Mass + SphP[i].StarMassFeed)) / SOLAR_METALLICITY;
#else
      double Zsol = 0;
#endif

      double Density = (P[i].Mass + SphP[i].StarMassFeed) / SphP[i].Volume;
      
      SphP[i].OpacityScaling[CH_DUST] = Zsol * Density / Units;

      Units = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm;

      double n_H2 = SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / SphP[i].Volume / (2 * PROTONMASS / All.cf_UnitMass_in_g);

      SphP[i].OpacityScaling[CH_H2] = n_H2 / Units;

      double n_Ionizing[3] = {SphP[i].GrackleSpeciesConserved(GRACKLE_HI) / SphP[i].Volume / (PROTONMASS / All.cf_UnitMass_in_g), 
                              SphP[i].GrackleSpeciesConserved(GRACKLE_HeI) / SphP[i].Volume / (4 * PROTONMASS / All.cf_UnitMass_in_g), 
                              SphP[i].GrackleSpeciesConserved(GRACKLE_HeII) / SphP[i].Volume / (4 * PROTONMASS / All.cf_UnitMass_in_g)};

      for(int s = 0; s < 3; s++)
        SphP[i].OpacityScaling[CH_HI + s] = n_Ionizing[s] / Units;
    }
}

double dtau_IR(int i, double length)
{
  double kappa_rerad = 1.0;

  double Dtau_IR = All.IRDtauMomentumBoostCoeff * kappa_rerard * SphP[i].OpacityScaling[CH_DUST] * length;

  return Dtau_IR;
}

static double H2Tab_A[H2TAB_N]; /* A at table nodes */
static double H2Tab_dlogN; /* log10 spacing */
static double H2Tab_A_thinmin; /* A(NMIN) = SigmaH2 * NMIN */

/* Wolcott-Green et al. (2011) self-shielding function */
static inline double f_selfshield_H2(double N_H2)
{
  double x = N_H2 / 5.0e14;
  double sq = sqrt(1.0 + x);
  return 0.965 / pow(1.0 + x / H2_SHIELD_B5, 1.1)
       + 0.035 / sq * exp(-8.5e-4 * sq);
}

/* Build A(N) once at startup (trapezoid, log-spaced with linear-N areas,
   16 sub-steps per interval so table error << fit error) */
void init_h2shield(void)
{
  H2Tab_dlogN = (H2TAB_LOGNMAX - H2TAB_LOGNMIN) / (H2TAB_N - 1);
  H2Tab_A_thinmin = SigmaH2 * pow(10.0, H2TAB_LOGNMIN);

  /* Thin part below NMIN: f_sh=1 */
  double A = H2Tab_A_thinmin;
  H2Tab_A[0] = A;

  for(int i = 1; i < H2TAB_N; i++)
    {
      double N0 = pow(10.0, H2TAB_LOGNMIN + (i - 1) * H2Tab_dlogN);
      double N1 = pow(10.0, H2TAB_LOGNMIN + i * H2Tab_dlogN);

      const int nsub = 16;
      double dN = (N1 - N0) / nsub;
      for(int k = 0; k < nsub; k++)
        {
          double Na = N0 + k * dN;
          A += 0.5 * (f_selfshield_H2(Na) + f_selfshield_H2(Na + dN)) * dN * SigmaH2;
        }

      H2Tab_A[i] = A;
    }
}

/* A(N): thin analytic below NMIN, clamp above NMAX, linear-in-logN inside */
static inline double h2shield_A(double N_H2)
{
  if(N_H2 <= 0.0)
    return 0.0;

  double logN = log10(N_H2);

  /* f_sh = 1 exactly */
  if(logN <= H2TAB_LOGNMIN)
    return SigmaH2 * N_H2;

  /* Lines exhausted */
  if(logN >= H2TAB_LOGNMAX)
    return H2Tab_A[H2TAB_N - 1];

  double u = (logN - H2TAB_LOGNMIN) / H2Tab_dlogN;
  int j = (int)u;
  double f = u - j;

  return H2Tab_A[j] * (1.0 - f) + H2Tab_A[j + 1] * f;
}

/* Exact per-cell line optical depth for a cell adding dN_H2 to a ray
   that has already accumulated N_H2 */
double h2shield_dtau(double N_H2, double dN_H2)
{
  double dtau = h2shield_A(N_H2 + dN_H2) - h2shield_A(N_H2);
  return dtau > 0.0 ? dtau : 0.0;
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

  R[0][0] = 1.0 - 2.0 * (y*y + z*z);
  R[0][1] = 2.0 * (x*y - z*w);
  R[0][2] = 2.0 * (x*z + y*w);

  R[1][0] = 2.0 * (x*y + z*w);
  R[1][1] = 1.0 - 2.0 * (x*x + z*z);
  R[1][2] = 2.0 * (y*z - x*w);

  R[2][0] = 2.0 * (x*z - y*w);
  R[2][1] = 2.0 * (y*z + x*w);
  R[2][2] = 1.0 - 2.0 * (x*x + y*y);
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
      cached_seed  = seed;
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
  free(w->rays); free(w);
}

/*
 * Rays are born inside their host cell, which the star density
 * machinery already identified through the SPH neighbour loop
 */
static void init_rays(RayWorkStack *work)
{
  double xtmp, ytmp, ztmp;

  double SQRT3 = sqrt(3.0);

  int ray_idx = 0;

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
            Radiated_Cell[w].Energy  += MechanicalFeedback->Radiated[w].Energy;
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
        unsigned long long rotation_seed = seed_rotation(P[host].ID);

        /* Loop over rays for this host */
        for(int iray = 0; iray < NRays; iray++)
          {
            RayPacket ray = {0};

            ray.star_id = P[host].ID; 

            ray.cell = host;

            ray.pos[0] = 0.0;
            ray.pos[1] = 0.0;
            ray.pos[2] = 0.0;

            ray.rotation_seed = rotation_seed;
            healpix_dir(rotation_seed, NSIDE_MIN, iray, ray.dir);

            ray.t = 0.0;
            ray.t_maximum = SQRT3 * All.BoxSize;

            ray.nside = NSIDE_MIN;
            ray.healpix_pixel = iray;

            ray.active_bands = NO_IR_ACTIVE;

            for(int w = 0; w < WAVEBANDS; w++)
              {
                ray.Radiated[w].Energy  = Radiated_Cell[w].Energy / NRays;
                ray.Radiated[w].Photons = Radiated_Cell[w].Photons / NRays;

                ray.Radiated_Init[w].Energy  = Radiated_Cell[w].Energy / NRays;
                ray.Radiated_Init[w].Photons = Radiated_Cell[w].Photons / NRays;

                if(ray.Radiated[w].Energy <= 0.0 && ray.Radiated[w].Photons <= 0.0)
                  ray.active_bands &= (uint8_t)(~(1u << w));
              }

            if(ray.active_bands == 0)
              continue;

            ray.N_H2 = 0.0;

            ray.ray_id = ray_idx++;
            ray.home_task = ThisTask;

            append_ray(work, &ray);
          }
      }

#else

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

        unsigned long long rotation_seed = seed_rotation(MechanicalFeedbackData->StarID);

        /* Loop over rays for this star */
        for(int iray = 0; iray < NRays; iray++)
          {
            RayPacket ray = {0};

            ray.star_id = MechanicalFeedbackData->StarID;

            ray.cell = host;

            ray.pos[0] = xrel[0];
            ray.pos[1] = xrel[1];
            ray.pos[2] = xrel[2];

            ray.rotation_seed = rotation_seed;
            healpix_dir(rotation_seed, NSIDE_MIN, iray, ray.dir);

            ray.t = 0.0;
            ray.t_maximum = SQRT3 * All.BoxSize;

            ray.nside = NSIDE_MIN;
            ray.healpix_pixel = iray;

            ray.active_bands = NO_IR_ACTIVE;

            for(int w = 0; w < WAVEBANDS; w++)
              {
                ray.Radiated[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
                ray.Radiated[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;

                ray.Radiated_Init[w].Energy = MechanicalFeedback->Radiated[w].Energy / NRays;
                ray.Radiated_Init[w].Photons = MechanicalFeedback->Radiated[w].Photons / NRays;

                if(ray.Radiated[w].Energy <= 0.0 && ray.Radiated[w].Photons <= 0.0)
                  ray.active_bands &= (uint8_t)(~(1u << w));
              }

            if(ray.active_bands == 0)
              continue;

            ray.N_H2 = 0.0;

            ray.ray_id = ray_idx++;
            ray.home_task = ThisTask;

            append_ray(work, &ray);
          }
      }

#endif
    
    ev += SphP[host].Host;
  }
}

/* Splits to 4 child rays. Children inherit position, cell and path length,
   so they simply restart the exit search in the cell the parent entered. */
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

      healpix_dir(parent->rotation_seed, new_nside, children[k].healpix_pixel, children[k].dir);

      for(int w = 0; w < WAVEBANDS; w++)
        {
          children[k].Radiated[w].Energy = parent->Radiated[w].Energy * 0.25;
          children[k].Radiated[w].Photons = parent->Radiated[w].Photons * 0.25;

          children[k].Radiated_Init[w].Energy = parent->Radiated_Init[w].Energy * 0.25;
          children[k].Radiated_Init[w].Photons = parent->Radiated_Init[w].Photons * 0.25;
        }
    }
}

static RayExportBuffer *init_export_buffer(long long capacity)
{
  RayExportBuffer *buf = malloc(sizeof(RayExportBuffer));

  buf->n = 0;
  buf->capacity = capacity;
  buf->task = malloc(capacity * sizeof(int));
  buf->rays = malloc(capacity * sizeof(RayPacket));

  return buf;
}

void append_export(RayExportBuffer *buf, const RayPacket *ray, int task)
{
  if(buf->n >= buf->capacity)
    {
      buf->capacity *= 2;
      buf->task = realloc(buf->task, buf->capacity * sizeof(int));
      buf->rays = realloc(buf->rays, buf->capacity * sizeof(RayPacket));
    }

  buf->task[buf->n] = task;
  buf->rays[buf->n] = *ray;
  buf->n++;
}

static void free_export_buffer(RayExportBuffer *buf)
{
  free(buf->rays); free(buf->task); free(buf);
}

static void sort_by_task(RayExportBuffer *buf)
{
  int *count = calloc(NTask, sizeof(int));
  int *cursor = malloc(NTask * sizeof(int));

  for(long long i = 0; i < buf->n; i++)
    count[buf->task[i]]++;

  cursor[0] = 0;
  for(int t = 1; t < NTask; t++)
    cursor[t] = cursor[t-1] + count[t-1];

  int *sorted_task = malloc(buf->n * sizeof(int));
  RayPacket *sorted_rays = malloc(buf->n * sizeof(RayPacket));

  for(long long i = 0; i < buf->n; i++)
    {
      int t = buf->task[i];
      long long pos = cursor[t]++;
      sorted_task[pos] = t;
      sorted_rays[pos] = buf->rays[i];
    }

  memcpy(buf->task, sorted_task, buf->n * sizeof(int));
  memcpy(buf->rays, sorted_rays, buf->n * sizeof(RayPacket));

  free(sorted_rays);
  free(sorted_task);
  free(cursor);
  free(count);
}

static void exchange_rays(RayExportBuffer *send, RayWorkStack *work)
{
  static MPI_Datatype MPI_RAYPACKET = MPI_DATATYPE_NULL;
  if(MPI_RAYPACKET == MPI_DATATYPE_NULL)
    {
      MPI_Type_contiguous(sizeof(RayPacket), MPI_BYTE, &MPI_RAYPACKET);
      MPI_Type_commit(&MPI_RAYPACKET);
    }

  /* Heap, not VLA: int send_count[NTask] blows the C stack for large NTask */
  int *send_count = malloc(NTask * sizeof(int));
  int *recv_count = malloc(NTask * sizeof(int));
  int *send_offset = malloc(NTask * sizeof(int));
  int *recv_offset = malloc(NTask * sizeof(int));

  memset(send_count, 0, NTask * sizeof(int));
  for(long long i = 0; i < send->n; i++)
    send_count[send->task[i]]++;

  MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  send_offset[0] = recv_offset[0] = 0;
  long long total_recv = recv_count[0];
  for(int i = 1; i < NTask; i++)
    {
      send_offset[i] = send_offset[i-1] + send_count[i-1];
      recv_offset[i] = recv_offset[i-1] + recv_count[i-1];
      total_recv += recv_count[i];
    }

  while(work->n + total_recv > work->capacity)
    {
      work->capacity *= 2;
      work->rays = realloc(work->rays, work->capacity * sizeof(RayPacket));
    }

  /* In ray Units */
  sort_by_task(send);

  MPI_Alltoallv(send->rays, send_count, send_offset, MPI_RAYPACKET,
  work->rays + work->n, recv_count, recv_offset, MPI_RAYPACKET,
  MPI_COMM_WORLD);

  work->n += total_recv;

  free(recv_offset);
  free(send_offset);
  free(recv_count);
  free(send_count);
}

static void radiation_feedback(void)
{
  /* Indexed by ionizing species a = 0,1,2 (HI, HeI, HeII) */
  static const double IonThreshold_eV[3] = {13.6, 24.6, 54.4};
  static const int IonGrackle[3] = {GRACKLE_HI, GRACKLE_HeI, GRACKLE_HeII};
  static const double IonAtomicMass[3] = {1.0, 4.0, 4.0};

  const double L3 = All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm * All.cf_UnitLength_in_cm;

  for(int i = 0; i < NumGas; i++)
    {
      if(P[i].Type != 0 || P[i].Mass == 0 || P[i].ID == 0)
        continue;

      const double V = SphP[i].Volume;
      const double dt = (P[i].TimeBinHydro ? (((integertime)1) << P[i].TimeBinHydro) : 0) * All.Timebase_interval;

      if(dt <= 0.0 || V <= 0.0)
        goto reset;

      const double V_cgs  = V * L3;
      const double dt_cgs = dt * All.cf_UnitTime_in_s;

#ifdef PHOTOELECTRIC_HEATING
      /* TrueAbsorbedFraction already applied in ray_deposit */
      const double epsilon_pe = 0.05;
      const double E_pe = SphP[i].AbsorbedPE * epsilon_pe * All.cf_UnitEnergy_in_cgs;

      SphP[i].PE_VolHeatingRate += E_pe / dt_cgs / V_cgs;
#endif

#ifdef DISSOCIATION
      const double n_H2 = SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) / V
                          / (2.0 * PROTONMASS / All.cf_UnitMass_in_g);

      /* AbsorbedH2Line holds pumped photons; F_DISS is the branching ratio */
      if(n_H2 > 0.0)
        SphP[i].H2_DissociationRate += (F_DISS * SphP[i].AbsorbedH2Line / (dt / All.cf_hubble_a / All.HubbleParam) / V) / n_H2;
#endif

#ifdef PHOTOIONIZATION
      for(int s = 0; s < 3; s++)
        {
          const double n = SphP[i].GrackleSpeciesConserved(IonGrackle[s]) / V
                           / (IonAtomicMass[s] * PROTONMASS / All.cf_UnitMass_in_g);

          const double N_abs = SphP[i].AbsorbedIonizing[s].Photons;
          const double E_abs = SphP[i].AbsorbedIonizing[s].Energy * All.cf_UnitEnergy_in_cgs;

          const double E_exc = E_abs - N_abs * IonThreshold_eV[s] * ELECTRONVOLT_IN_ERGS;

          if(n <= 0.0)
            continue;

          const double n_cgs = n / L3;

          if(E_exc > 0.0)
            SphP[i].IonHeatingRate[s] += (E_exc / dt_cgs / V_cgs) / n_cgs;
          else if(N_abs > 0.0)
            warn("STAR_RADIATION: sub-threshold mean photon energy, species %d, cell %d "
                 "(E_abs=%g N_abs=%g) \n", s, i, E_abs, N_abs);

          SphP[i].IonizationRate[s] += (N_abs / (dt / All.cf_hubble_a / All.HubbleParam) / V) / n;
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
      for(int s = 0; s < 3; s++)
        SphP[i].AbsorbedIonizing[s].Energy = SphP[i].AbsorbedIonizing[s].Photons = 0.0;
#endif
    }
}

#ifdef PHOTOIONIZATION
static void rt_timestep(void)
{
  int idx, i;
  double eps_ion = All.RTIonizationTimestepFraction;

  for(idx = 0; idx < TimeBinsHydro.NActiveParticles; idx++)
    {
      i = TimeBinsHydro.ActiveParticleList[idx];
      if(i < 0)
        continue;

      /* Normalize ionization rate by total hydrogen (was normalized by neutral for grackle) */
      double m_HI = SphP[i].GrackleSpeciesConserved(GRACKLE_HI);

      double m_H = SphP[i].GrackleSpeciesConserved(GRACKLE_HI) + SphP[i].GrackleSpeciesConserved(GRACKLE_HII);

#if GRACKLE_CHEMISTRY >= 2
      m_H += SphP[i].GrackleSpeciesConserved(GRACKLE_H2I) + SphP[i].GrackleSpeciesConserved(GRACKLE_H2II)
          + SphP[i].GrackleSpeciesConserved(GRACKLE_HM);
#endif

#if GRACKLE_CHEMISTRY >= 3
      m_H += SphP[i].GrackleSpeciesConserved(GRACKLE_DI) + SphP[i].GrackleSpeciesConserved(GRACKLE_DII)
          + SphP[i].GrackleSpeciesConserved(GRACKLE_HDI);
#endif

      double rate = 0.0;
      if(m_H > 0)
        {
          double x_HI = m_HI / m_H;
          rate = SphP[i].IonizationRate[0] * x_HI;
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
      for(int s = 0; s < 3; s++)
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
  long long export_capacity = n_rays_local > 0 ? n_rays_local : 1024;

  RayWorkStack *work = init_work_stack(work_capacity);
  RayExportBuffer *export_buf = init_export_buffer(export_capacity);

  init_rays(work);

  long long n_global;
  int iter = 0;
  do
    {
      t0 = second();

      /* Do local work first then exchange */
      while(work->n > 0)
        {
          RayPacket ray = work->rays[--work->n];
          raytrace_voronoi(&ray, work, export_buf);
        }

      /* Send rays to remote ranks, receive rays from remote ranks */
      exchange_rays(export_buf, work);

      /* Reset export buffer for this round */
      export_buf->n = 0;

      /* Check if anyone still has rays in flight */
      sumup_longs(1, &work->n, &n_global);

      t1 = second();

      iter++;

      if(n_global > 0 && iter > 0)
        mpi_printf("STAR_RADIATION: Rad iteration %3d: need to repeat for %12lld rays. (took %g sec)\n", iter, n_global,
        timediff(t0, t1));

      if(iter > 4 * MAXITER)
        terminate("STAR_RADIATION: %lld rays still in flight after %d iterations\n", n_global, iter);

    } while(n_global > 0);

  radiation_feedback();

#ifdef PHOTOIONIZATION
  rt_timestep();
#endif

  free_export_buffer(export_buf);
  free_work_stack(work);

  TIMER_STOP(CPU_STARS_RADIATION);
}