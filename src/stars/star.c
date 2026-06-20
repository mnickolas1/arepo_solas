#include "../main/allvars.h"


int NumStars;

#ifdef STAR_FEEDBACK_ACTIVE
struct TimeBinData TimeBinsStar;
#endif

#if defined(WINDS) || defined(SUPERNOVAE)
Mechanical_Feedback_Pack MechanicalFeedbackEvents;
#endif

Star_Particle_Data *SP;
