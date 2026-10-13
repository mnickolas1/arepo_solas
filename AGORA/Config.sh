#---------------------------------------- SOLAS additionsS

#---------------------------------------- Cooling parameters
COOLING                # General cooling flag

USE_GRACKLE            # Enable grackle
GRACKLE_CHEMISTRY=2    # Set grackle chemistry mode 

#---------------------------------------- Star options
STARS                  # General stars framework flag

STAR_PARTICLES=2       # Star particles model flag: set to 0, 1 for massive star particles, set to 2 for resolved individual stars

STAR_IN_CELL           # Place stars on the host cell generator; superpose all stars in a cell into one source

PHOTOIONIZATION        # Only photoionization

#---------------------------------------- Radiation / RT options
RAD_TOTAL_TRUNCATION   # Truncate rays against the initial total energy/photons, not the per-band values
RT_TIMESTEP            # RT timestep limiter (requires DISSOCIATION or PHOTOIONIZATION)

RT_STATISTICS          # Ray-tracing statistics
RT_COMM_STATISTICS     # RT communication timing statistics

DOMAIN_MEASURED_WEIGHTS # Use measured per-channel timings for the domain decomposition weights

RTYPE
#DTYPE

#RT_COMM_SYNC

#--------------------------------------- Arepo public
CREATE_FULL_MESH

FORCE_EQUAL_TIMESTEPS

#--------------------------------------- Mesh motion and regularization; default: moving mesh
REGULARIZE_MESH_CM_DRIFT      # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE    # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Refinement and derefinement; default: no refinement/derefinement; criterion: target mass
REFINEMENT_SPLIT_CELLS        # Refinement
REFINEMENT_MERGE_CELLS        # Derefinement

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS          # non-local timestep criterion (take 'signal speed' into account)

#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1             # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.

#--------------------------------------- output options
HAVE_HDF5                     # needed when HDF5 I/O support is desired (recommended)

#--------------------------------------- Testing and Debugging options
DEBUG                         # enables core-dumps
