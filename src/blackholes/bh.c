#include "../stars/star.h"

int NumBhs;

FILE *FdBlackHoles; 

#if defined(BH_ACCRETION_ACTIVE) || defined(BH_FEEDBACK_ACTIVE)
struct TimeBinData TimeBinsBh;
#endif

struct bh_particle_data *BhP;