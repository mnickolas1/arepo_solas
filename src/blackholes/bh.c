#include "../main/allvars.h"


int NumBhs;

FILE *FdBlackHoles; 

#ifdef BH_ACTIVE
struct TimeBinData TimeBinsBh;
#endif

struct Bh_Particle_Data *BhP;