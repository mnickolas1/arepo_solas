#include "../stars/star.h"

int NumStars;

struct star_particle_data *SP;

#ifdef STAR_FEEDBACK_ACTIVE
struct TimeBinData TimeBinsStar;
#endif