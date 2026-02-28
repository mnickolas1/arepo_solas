#define MAX_RAYS 12288 // enough for NSIDE=32, adjust as needed

extern double HealpixDirs[MAX_RAYS][3];
extern int NRays; // 12 * NSIDE^2

extern double RayColumnDensity[MAX_RAYS];
extern double RayIntensity[MAX_RAYS];