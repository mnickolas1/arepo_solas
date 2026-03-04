#define RAD_BACKGROUND 1 
#define MAX_RAYS 12288 
#define RAY_STACK_SIZE 64
#define MAX_RAYS_EXCHANGE 10000

extern double HealpixDirs[MAX_RAYS][3];
extern int NRays; // 12 * NSIDE^2

typedef struct 
{
  double t_enter;
  int node;
} StackEntry;

typedef struct 
{
  double pos[3];
  double dir[3];
  double RAD_Ionizing;
  int ray_id;
  int home_task;
  
  /* pending top-level nodes still to traverse after current domain */
  StackEntry pending[RAY_STACK_SIZE];
  int n_pending;
  int target_node;
  double t;
} RayData;

typedef struct 
{
  int n; /* number of rays to export */
  RayData *rays; /* ray data */
  int *task; /* which task to send each ray to */
  int capacity; /* allocated capacity */
} RayExportBuffer;

extern struct rad_resultsactiveimported_data
{
  double RAD_Ionizing;
  int index; /* local SphP index on home task */
} *Rad_ResultsActiveImported;

/* Wavebands */
typedef enum
{
  LYMAN_WERNER = 0,
  ULTRAVIOLET,
  OPTICAL,
  INFRARED,
  WB_COUNT
} Waveband;

/* Opacity coefficients */
extern double Kappa[WB_COUNT] = {
  1.0,  // LW
  1.0,  // UV
  1.0,  // OP
  1.0   // IR
};