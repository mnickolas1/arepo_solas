/* star functions */

/* Memory allocation */
void reallocate_memory_maxpartstars(void);
void domain_resize_storage_stars(int count_get_star);

/* Timesteps */
integertime get_timestep_star(int p);
void update_star_timesteps(void);
void reconstruct_star_timebins(void);
void update_list_of_active_star_particles(void);

/* Density-Feedback loop */
void star_density(void);
void star_prep(void);
void star_feedback(void);
void perform_end_of_step_star_physics(void);
void star_kernel(double u, double hinv3, double hinv4, double *wk, double *dwk);

/* Feedback tables interpolation */
void load_stellar_tables(const char *filename);
void free_stellar_tables(void);
struct star_feedback star_feedback_compute(double dt, double z_val, double m_val, double a);
struct star_feedback units_for_feedback(struct star_feedback star);

double IntegralTrapezoidal(double a, double b, int N, double (*f)(double));

#ifdef STAR_PARTICLES
/* IMF */
double imf_kroupa(double m); 
double imf_chabrier(double m); 
double imf_salpeter(double m);
double imf(double m); 
double m_times_imf(double m); 
void build_imf_cdf(void);
double sample_imf(double u);
#endif

#if STAR_PARTICLES = 1
void setup_mass_bins(void);
void sample_star_particle(double m, int *bins);
struct star_feedback star_particle_feedback(int index, double dt, double z, double a);
#endif

#ifdef INDIVIDUAL_STAR_BY_STAR_FORMATION
void sf_starbystar(void);
void sf_massdrain(void);
#endif