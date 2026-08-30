#include "../main/allvars.h"
#include "../main/proto.h"


int NumStars;

#ifdef STAR_FEEDBACK_ACTIVE
struct TimeBinData TimeBinsStar;
#endif

#ifdef STAR_FEEDBACK_ACTIVE
Mechanical_Feedback_Events MechanicalFeedbackEvents;
#endif

Star_Particle_Data *SP;

void star_init(void)
{
#ifdef STAR_PARTICLES
  if(ThisTask == 0)
    {
      build_imf_cdf();

#if STAR_PARTICLES < 2
      setup_mass_bins();
#endif

#if STAR_PARTICLES == 0
      setup_imf_integrals();
#endif
    }
  
  MPI_Bcast(cdf_masses, N_CDF_BINS + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(cdf_values, N_CDF_BINS + 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#if STAR_PARTICLES < 2
  MPI_Bcast(StarMeanMassInBins, NBINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#if STAR_PARTICLES == 0
  MPI_Bcast(&norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(bin_imf, NBINS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#if STAR_PARTICLES < 2
  if(RestartFlag != 1)
    {
      for(int i = 0; i < NumStars; i++)
        sample_star_particle(PPS(i).Mass * All.cf_UnitMass_in_Msun, SP[i].NumOfStarsInBins);
    }
#endif
#endif

#ifdef STAR_FEEDBACK_ACTIVE  
  mpi_printf("Loading star evolution tables\n");
  load_star_tables(All.StarTablesFile);

  feedback_init(&MechanicalFeedbackEvents);

  if(RestartFlag != 1)
    {
      for(int i = 0; i < NumStars; i++)
        SP[i].WithFeedback = 1;
    }
#endif

#ifdef STAR_RADIATION_ACTIVE
  init_h2shield();
#endif
}