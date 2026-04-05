#ifndef BH_PROTO_H
#define BH_PROTO_H

/* black hole functions */

/* Memory allocation */
void reallocate_memory_maxpartbhs(void);
void domain_resize_storage_bhs(int count_get_bh);

#if defined(BH_ACCRETION_ACTIVE) || defined(BH_FEEDBACK_ACTIVE)
/* Timesteps */
integertime get_timestep_bh(int p);
void update_bh_timesteps(void);
void reconstruct_bh_timebins(void);
void update_list_of_active_bh_particles(void);

/* Density loop */
void bh_density(void);
#endif

#ifdef BH_ACCRETION_ACTIVE
/* Accretion loops */
void bh_accretion(void);
void bh_swallow(void);
#endif

/* Feedback loops */
#ifdef BH_THERMAL_FEEDBACK
void bh_feedback(void);
#endif

#ifdef BH_JET_FEEDBACK
void bh_jet_density(void);
void bh_jet_feedback(void);
#endif

void blackhole_mark_cells_for_refinement(void);

#endif 