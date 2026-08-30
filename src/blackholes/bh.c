#include "../main/allvars.h"


int NumBhs;

FILE *FdBlackHoles; 

#ifdef BH_ACTIVE
struct TimeBinData TimeBinsBh;
#endif

struct Bh_Particle_Data *BhP;

void bh_init(void);
{
#ifdef BH_FEEDBACK_ACTIVE
  if(RestartFlag != 1) 
    {
      All.FeedbackFlag = 1;
    }
#endif
} 