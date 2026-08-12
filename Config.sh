#--------------------------------------- SOLAS additions

#--------------------------------------- Cooling parameters
USE_GRACKLE
GRACKLE_CHEMISTRY=1    

#--------------------------------------- Star options
STARS                  # General stars framework flag
STAR_PARTICLES=2       # Star particles model flag: set to 0, 1 for massive star particles, set to 2 for resolved individual stars

PHOTOIONIZATION      


#--------------------------------------- Arepo public
# rtype
NOHYDRO

CREATE_FULL_MESH

FORCE_EQUAL_TIMESTEPS


#--------------------------------------- Mesh motion and regularization; default: moving mesh
REGULARIZE_MESH_CM_DRIFT      # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE    # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Refinement and derefinement; default: no refinement/derefinement; criterion: target mass
REFINEMENT_SPLIT_CELLS        # Refinement
REFINEMENT_MERGE_CELLS        # Derefinement

#--------------------------------------- non-standard phyiscs
COOLING                       # Simple primordial cooling (need flag to call grackle -> overrides)

#--------------------------------------- Gravity treatment; default: no gravity
GRAVITY_NOT_PERIODIC          # gravity is not treated periodically

#--------------------------------------- Single/Double Precision
DOUBLEPRECISION=1              # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION
#NGB_TREE_DOUBLEPRECISION      # if this is enabled, double precision is used for the neighbor node extension

#--------------------------------------- output options
HAVE_HDF5                     # needed when HDF5 I/O support is desired (recommended)

#--------------------------------------- Testing and Debugging options
DEBUG                         # enables core-dumps

OVERRIDE_PEANOGRID_WARNING