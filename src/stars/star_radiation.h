#ifndef STAR_RADIATION_H
#define STAR_RADIATION_H

#include <stdint.h>


/*
 *   RADIATION_PRESSURE      direct + IR-reradiated momentum coupling
 *   PHOTOELECTRIC_HEATING   grain photoelectric heating (FUV: UV + LW dust)
 *   DISSOCIATION            H2 Lyman-Werner photodissociation (+ self-shielding)
 *   PHOTOIONIZATION         HI / HeI / HeII photoionization + photoheating
 */

#define RAD_TRUNC_FRAC 1.0e-2

/* Below this optical depth -expm1(-tau) is evaluated as a series */
#define RAD_TAU_THIN 1.0e-3

/* Safety cap on the number of Voronoi cells a single ray may cross in one call to raytrace_voronoi() */ 
/* Only trips on degenerate geometry */
#define RAY_MAX_CELL_STEPS (1 << 22)

/* Relative tolerance (in units of the host cell radius) used to break ties
   between coincident bisectors at cell edges/vertices */
#define RAY_TOL 1.0e-10

/* Healpix: Powers of 2! */
#define NSIDE_MIN 1
#define NSIDE_MAX 128

/* Starting number of rays */
#define NRays (12 * NSIDE_MIN * NSIDE_MIN)

/* Dissociation of H2 */
#define SIGMA_DISS 2.47e-18 /* cm^2, dissociation-weighted eff. cross-section (DB96, Baczynski+15) */

#define F_DISS 0.15 /* dissociation branching per absorption */

#define H2_SHIELD_B5 3.0 /* Doppler b / (km/s); fixed */

/* Shielding log table parameters */
#define H2TAB_N 1024
#define H2TAB_LOGNMIN 11.0 /* log10 N_H2 [cm^-2]: f_sh = 1 below */
#define H2TAB_LOGNMAX 24.0 /* f_sh negligible above */

/*
 * INFRARED: inf A - 12398.4 A (0 eV - 1 eV)
 *
 * OPTICAL: 12398.4 A - 2066.4 A (1 eV - 6 eV)
 *
 * ULTRAVIOLET: 2066.4 A - 1107.0 A (6 eV - 11.2 eV)
 *
 * LYMAN_WERNER: 1107.0 A - 911.6 A (11.2 eV - 13.6 eV)
 *
 * IONIZING_HI: 911.6 A - 504.0 A (13.6 eV - 24.6 eV)
 *
 * IONIZING_HeI: 504.0 A - 227.9 A (24.6 eV - 54.4 eV)
 *
 * IONIZING_HeII: 227.9 A - 0 A (54.4 eV - inf eV)
 */
typedef enum
{ 
  INFRARED = 0,
  OPTICAL,
  ULTRAVIOLET,
  LYMAN_WERNER,
  IONIZING_HI,
  IONIZING_HeI,
  IONIZING_HeII,
  WAVEBANDS
} Waveband;

enum
{
  CH_DUST = 0,
  CH_H2,
  CH_HI,
  CH_HeI,
  CH_HeII,
  CHANNELS
};

/* Which channels can absorb in each band */
static const uint8_t BandChannels[WAVEBANDS] =
{
  [INFRARED] = (1u << CH_DUST),
  [OPTICAL] = (1u << CH_DUST),
  [ULTRAVIOLET] = (1u << CH_DUST),
  [LYMAN_WERNER] = (1u << CH_DUST) | (1u << CH_H2),
  [IONIZING_HI] = (1u << CH_DUST) | (1u << CH_HI),
  [IONIZING_HeI] = (1u << CH_DUST) | (1u << CH_HI) | (1u << CH_HeI),
  [IONIZING_HeII] = (1u << CH_DUST) | (1u << CH_HI) | (1u << CH_HeI) | (1u << CH_HeII),
};

/* Active bands - change at init_rays */
#define ALL_BANDS_ACTIVE ((uint8_t)((1u << WAVEBANDS) - 1u))
#define NO_IR_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ~(1u << INFRARED)))
#define NO_IONIZING_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ~(1u << IONIZING_HI) & ~(1u << IONIZING_HeI) & ~(1u << IONIZING_HeII)))
#define ONLY_IONIZING_ACTIVE ((uint8_t)(ALL_BANDS_ACTIVE & ((1u << IONIZING_HI) | (1u << IONIZING_HeI) | (1u << IONIZING_HeII))))

/* Bands carrying photons */
static const uint8_t BandTrackPhotons = (1u << LYMAN_WERNER) | (1u << IONIZING_HI) | (1u << IONIZING_HeI) | (1u << IONIZING_HeII);

typedef struct WavebandData
{
  double Energy;
  double Photons;
} WavebandData;

typedef struct
{
  double E[CHANNELS];
  double N[CHANNELS];
} ChannelsDtau;

typedef struct
{
  uint8_t mask;
  WavebandData Band[WAVEBANDS]; /* total removed from the ray */
  WavebandData Ch[WAVEBANDS][CHANNELS]; /* attribution; sums to Band */
} Absorption;

/* Dust */
extern double Kappa_E[WAVEBANDS];
extern double Kappa_N[WAVEBANDS];

extern double TrueAbsorbedFraction[WAVEBANDS];
extern double ReradiatedFraction[WAVEBANDS];

/* H2 lines */
extern double SigmaH2;

/* Ionizing */
extern double Sigma_E[3][3];
extern double Sigma_N[3][3];

/*
 * RayPacket
 * ---------
 * A ray is fully described by the cell it currently occupies plus its
 * entry point into that cell. Because the Voronoi cell is exactly the
 * intersection of the half spaces
 *
 *     (x - m_ij) . d_ij <= 0 ,   d_ij = s_j - s_i ,  m_ij = (s_i + s_j)/2
 *
 * over the face-defining neighbours listed in DC(i), the exit point and the
 * next cell follow from a single min-reduction over that list.
 */
typedef struct RayPacket
{
  MyIDType star_id;
    
  /* Cell currently occupied: local SphP index on the owning task */
  int cell;

  double pos[3]; /* Position relative to the generator of `cell` */
  double dir[3]; /* Unit propagation direction (global frame) */

  double t; /* Path length accumulated since the star */
  double t_maximum; /* Hard stop */

  int nside; /* Current HEALPix nside level */
  int healpix_pixel; /* Current HEALPix pixel (NESTED) */
  uint8_t locate_head; /* Locate flag for split rays */

  /* Bitmask: bit w is SET while band w is still alive */
  /* When active_bands == 0 the ray is fully absorbed - return immediately */
  uint8_t active_bands;

  /* Ray energy and photons */
  WavebandData Radiated[WAVEBANDS];

#ifdef RAD_TOTAL_TRUNCATION
  /* Sum over bands at birth */
  double E_init; 
  double N_init; 
#else
  WavebandData Radiated_Init[WAVEBANDS];
#endif

  /* Accumulated H2 column since source */
  double N_H2; /* cgs! */
} RayPacket;

typedef struct RayWorkStack
{
  long long n; /* Number of rays on this stack */
  long long capacity; /* Allocated capacity */
  RayPacket *rays;
} RayWorkStack;

typedef struct RayExportBuffer
{
  long long n; /* Number of rays to export */
  long long capacity; /* Allocated capacity */
  int *ngbs; /* Ngbs slot */
  RayPacket *rays;
} RayExportBuffer;

extern int RayNgbNTask; /* number of mesh-neighbour ranks */
extern int *RayNgbTask; /* ascending list of neighbour ranks, length RayNgbNTask */
extern int *RayTaskToNgb; /* rank -> neighbour slot, or -1; length NTask */

/* Cell steps between in-traversal progress calls */
/* A single ray may cross millions of cells; without this a long traversal stalls */
#define RAY_STEPS_PROGRESS_MASK 1023

#ifdef RT_COMM_SYNC
 
typedef RayExportBuffer RayComms;
 
/* No-op: the synchronous path only communicates between rounds */
static inline void ray_comms_progress(RayComms *comm) { (void)comm; }
 
#else
 
/* Packets per message */
#define RAY_MSG_MAX 256
 
/* Send slots beyond the one-per-neighbour minimum */
#define RAY_SEND_SPARE 32
 
/* Pre-posted receives: 2 per neighbour */
#define RAY_RECV_SLOTS_MIN 8
#define RAY_RECV_SLOTS_MAX 64
 
/* Rays traced between calls into the comm layer */
#define RAY_PROGRESS_CHUNK 64
 
/* Progress calls between forced flushes of partially filled send buffers */
#define RAY_FLUSH_INTERVAL 8
 
/* Spin count at which the driver starts complaining about a possible hang */
#define RAY_ASYNC_SPIN_WARN 50000000LL
 
typedef struct RayCommsAsync RayComms;
#endif

#endif