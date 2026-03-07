#ifndef STAR_RADIATION_H
#define STAR_RADIATION_H

#include <stdint.h>

#define RAD_BACKGROUND 0
#define RAD_TRUNC_FRAC 0.01 
#define MAX_NUM_RAYS 12288 
#define RAY_STACK_SIZE 64
#define MAX_NUM_RAYS_TO_EXCHANGE 10000

#define ALL_BANDS_ACTIVE  ((uint8_t)((1u << WAVEBANDS) - 1u))

extern double HealpixDirs[MAX_NUM_RAYS][3];
extern int NRays; // 12 * NSIDE^2

typedef enum
{ IONIZING = 0,
  LYMAN_WERNER,
  ULTRAVIOLET,
  OPTICAL,
  INFRARED,
  WAVEBANDS
} Waveband;

extern double Kappa[WAVEBANDS];

typedef struct 
{
  double t_enter;
  int node;
} StackEntry;

typedef struct 
{
  double pos[3];
  double dir[3];
  double RAD[WAVEBANDS];
  double RAD_Initial[WAVEBANDS];

  /* Bitmask: bit w is SET while band w is still alive.
     Cleared when RAD[w] < RAD_TRUNC_FRAC * RAD_Initial[w].
     When active_bands == 0 the ray is fully absorbed – return immediately. */
  uint8_t  active_bands;

  int ray_id;
  int home_task;
  
  /* pending top-level nodes still to traverse after current domain */
  StackEntry pending[RAY_STACK_SIZE + 1];
  int n_pending;
  int target_node;
  double t;
} RayPacket;

typedef struct 
{
  int n; /* number of rays to export */
  RayPacket *rays; /* ray information */
  int *task; /* which task to send each ray to */
  int capacity; /* allocated capacity */
} RayExportBuffer;

extern struct rad_resultsactiveimported_data
{
  double RAD[WAVEBANDS];
  int index; /* local SphP index on home task */
} *Rad_ResultsActiveImported;

#endif