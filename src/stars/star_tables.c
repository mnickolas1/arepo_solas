#include <hdf5.h>

#include "../main/allvars.h"
#include "../main/proto.h"


/* Main Sequence */
int Z_COUNT = 0;
int M_COUNT = 0;

double *Z_VALUES = NULL;
double *M_VALUES = NULL;

double *logZ_VALUES = NULL;
double *logM_VALUES = NULL;
   
int **N = NULL;

double ***Age = NULL;
double ***FractionalAge = NULL;
double ***logRadius = NULL;
double ***logTemperature = NULL;

#ifdef WINDS
/* Winds */
double ***logMassLossRate = NULL;

#if GRACKLE_CHEMISTRY >= 1
double ***WindX = NULL;
double ***WindY = NULL;
#endif
#ifdef METALS
double ***WindZ = NULL;
#endif
double ***logWindVelocity = NULL;
#endif

#ifdef STAR_RADIATION_ACTIVE
/* Radiation */
WavebandData ***logFlux[WAVEBANDS] = {0};
#endif 

#ifdef SUPERNOVAE
/* Supernovae */ 
double **SN_MassLoss = NULL;
#if GRACKLE_CHEMISTRY >= 1
double **SN_X = NULL;
double **SN_Y = NULL;
#endif  
#ifdef METALS  
double **SN_Z = NULL;
#endif
#endif  

#ifdef AGB 
/* Asymptotic Giant Branch */
double **AGB_MassLoss; 
#ifdef METALS
double **AGB_MetalsLoss; 
#endif 
#endif

void free_stellar_tables(void)
{
  if(N)
    {
      for(int z = 0; z < Z_COUNT; z++)
        for(int m = 0; m < M_COUNT; m++)
          {
            free(Age[z][m]);
            free(FractionalAge[z][m]);
            free(logRadius[z][m]);
            free(logTemperature[z][m]);

#ifdef WINDS
            free(logMassLossRate[z][m]);
#if GRACKLE_CHEMISTRY >= 1
            free(WindX[z][m]);
            free(WindY[z][m]);
#endif
#ifdef METALS
            free(WindZ[z][m]);
#endif
            free(logWindVelocity[z][m]);
#endif

#ifdef STAR_RADIATION_ACTIVE
            for(int w = 0; w < WAVEBANDS; w++)
              free(logFlux[w][z][m]);
#endif
          }

      for(int z = 0; z < Z_COUNT; z++)
        {
          free(N[z]);
          
          free(Age[z]);
          free(FractionalAge[z]);
          free(logRadius[z]);
          free(logTemperature[z]);

#ifdef WINDS
          free(logMassLossRate[z]);
#if GRACKLE_CHEMISTRY >= 1
          free(WindX[z]);
          free(WindY[z]);
#endif
#ifdef METALS
          free(WindZ[z]);
#endif
          free(logWindVelocity[z]);
#endif

#ifdef STAR_RADIATION_ACTIVE
          for(int w = 0; w < WAVEBANDS; w++)
            free(logFlux[w][z]);
#endif

#ifdef SUPERNOVAE
          free(SN_MassLoss[z]);
#if GRACKLE_CHEMISTRY >= 1
          free(SN_X[z]);
          free(SN_Y[z]);
#endif
#ifdef METALS
          free(SN_Z[z]);
#endif
#endif
        }

      free(Z_VALUES);
      free(M_VALUES);

      free(logZ_VALUES);
      free(logM_VALUES);
      
      free(N);
      
      free(Age);
      free(FractionalAge);
      free(logRadius);
      free(logTemperature);

#ifdef WINDS
      free(logMassLossRate);
#if GRACKLE_CHEMISTRY >= 1
      free(WindX);
      free(WindY);
#endif
#ifdef METALS
      free(WindZ);
#endif
      free(logWindVelocity);
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        free(logFlux[w]);
#endif

#ifdef SUPERNOVAE
      free(SN_MassLoss);
#if GRACKLE_CHEMISTRY >= 1
      free(SN_X);
      free(SN_Y);
#endif
#ifdef METALS
      free(SN_Z);
#endif
#endif

      Z_VALUES = NULL;
      M_VALUES = NULL;

      logZ_VALUES = NULL;
      logM_VALUES = NULL;
      
      N = NULL;
      
      Age = NULL;
      FractionalAge = NULL;
      logRadius = NULL;
      logTemperature = NULL;

#ifdef WINDS
      logMassLossRate = NULL;
#if GRACKLE_CHEMISTRY >= 1
      WindX = NULL;
      WindY = NULL;
#endif
#ifdef METALS
      WindZ = NULL;
#endif
      logWindVelocity = NULL;
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        logFlux[w] = NULL;
#endif    

#ifdef SUPERNOVAE
      SN_MassLoss = NULL;
#if GRACKLE_CHEMISTRY >= 1
      SN_X = NULL;
      SN_Y = NULL;
#endif
#ifdef METALS
      SN_Z = NULL;
#endif
#endif
    }
}

void load_star_tables(const char *filename)
{ 
  hid_t file_id = -1;

  if(ThisTask == 0)
    {
      file_id = my_H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);  
      
      hid_t zc = my_H5Aopen_name(file_id, "Z_COUNT");
      hid_t mc = my_H5Aopen_name(file_id, "M_COUNT");
      
      my_H5Aread(zc, H5T_NATIVE_INT, &Z_COUNT, "Z_COUNT", 1);
      my_H5Aread(mc, H5T_NATIVE_INT, &M_COUNT, "M_COUNT", 1);

      my_H5Aclose(zc, "Z_COUNT");
      my_H5Aclose(mc, "M_COUNT");

      Z_VALUES = malloc(Z_COUNT * sizeof(double));
      M_VALUES = malloc(M_COUNT * sizeof(double));

      logZ_VALUES = malloc(Z_COUNT * sizeof(double));
      logM_VALUES = malloc(M_COUNT * sizeof(double));
      
      hid_t zv = my_H5Dopen(file_id, "LOGZ_VALUES");
      hid_t mv = my_H5Dopen(file_id, "LOGM_VALUES");

      my_H5Dread(zv, H5T_NATIVE_DOUBLE, 
              H5S_ALL, H5S_ALL, H5P_DEFAULT, logZ_VALUES, "LOGZ_VALUES");
      my_H5Dread(mv, H5T_NATIVE_DOUBLE, 
              H5S_ALL, H5S_ALL, H5P_DEFAULT, logM_VALUES, "LOGM_VALUES");

      for(int z = 0; z < Z_COUNT; z++)
        {
          if(z > 0 && logZ_VALUES[z] <= logZ_VALUES[z - 1])
            terminate("LOGZ_VALUES not strictly increasing at z=%d (%g <= %g)", z, logZ_VALUES[z], logZ_VALUES[z - 1]);
          Z_VALUES[z] = pow(10.0, logZ_VALUES[z]);
        }

      for(int m = 0; m < M_COUNT; m++)
        {
          if(m > 0 && logM_VALUES[m] <= logM_VALUES[m - 1])
            terminate("LOGM_VALUES not strictly increasing at m=%d", m);
          M_VALUES[m] = pow(10.0, logM_VALUES[m]);
        }

      my_H5Dclose(zv, "LOGZ_VALUES");
      my_H5Dclose(mv, "LOGM_VALUES");
    }

  MPI_Bcast(&Z_COUNT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&M_COUNT, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(ThisTask != 0)
    {
      Z_VALUES = malloc(Z_COUNT * sizeof(double));
      M_VALUES = malloc(M_COUNT * sizeof(double));

      logZ_VALUES = malloc(Z_COUNT * sizeof(double));
      logM_VALUES = malloc(M_COUNT * sizeof(double));
    }

  MPI_Bcast(Z_VALUES, Z_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(M_VALUES, M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Bcast(logZ_VALUES, Z_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(logM_VALUES, M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  N = malloc(Z_COUNT * sizeof(int*));

  Age = malloc(Z_COUNT * sizeof(double**));
  FractionalAge = malloc(Z_COUNT * sizeof(double**));
  logRadius = malloc(Z_COUNT * sizeof(double**));
  logTemperature = malloc(Z_COUNT * sizeof(double**));

#ifdef WINDS
  logMassLossRate = malloc(Z_COUNT * sizeof(double**));
#if GRACKLE_CHEMISTRY >= 1
  WindX = malloc(Z_COUNT * sizeof(double**));
  WindY = malloc(Z_COUNT * sizeof(double**));
#endif
#ifdef METALS
  WindZ = malloc(Z_COUNT * sizeof(double**));
#endif
  logWindVelocity = malloc(Z_COUNT * sizeof(double**));
#endif

#ifdef STAR_RADIATION_ACTIVE
  for(int w = 0; w < WAVEBANDS; w++)
    logFlux[w] = malloc(Z_COUNT * sizeof(WavebandData**));
#endif

#ifdef SUPERNOVAE
  SN_MassLoss = malloc(Z_COUNT * sizeof(double *));
#if GRACKLE_CHEMISTRY >= 1
  SN_X = malloc(Z_COUNT * sizeof(double *));
  SN_Y = malloc(Z_COUNT * sizeof(double *));
#endif
#ifdef METALS
  SN_Z = malloc(Z_COUNT * sizeof(double *));
#endif
#endif 

  for(int z = 0; z < Z_COUNT; z++)
    {
      N[z] = malloc(M_COUNT * sizeof(int));

      Age[z] = malloc(M_COUNT * sizeof(double*));
      FractionalAge[z] = malloc(M_COUNT * sizeof(double*));
      logRadius[z] = malloc(M_COUNT * sizeof(double*));
      logTemperature[z] = malloc(M_COUNT * sizeof(double*));

#ifdef WINDS
      logMassLossRate[z] = malloc(M_COUNT * sizeof(double*));
#if GRACKLE_CHEMISTRY >= 1
      WindX[z] = malloc(M_COUNT * sizeof(double*));
      WindY[z] = malloc(M_COUNT * sizeof(double*));
#endif
#ifdef METALS
      WindZ[z] = malloc(M_COUNT * sizeof(double*));
#endif
      logWindVelocity[z] = malloc(M_COUNT * sizeof(double*));
#endif

#ifdef STAR_RADIATION_ACTIVE
      for(int w = 0; w < WAVEBANDS; w++)
        logFlux[w][z] = malloc(M_COUNT * sizeof(WavebandData*));
#endif

#ifdef SUPERNOVAE
      SN_MassLoss[z] = malloc(M_COUNT * sizeof(double));
#if GRACKLE_CHEMISTRY >= 1
      SN_X[z] = malloc(M_COUNT * sizeof(double));
      SN_Y[z] = malloc(M_COUNT * sizeof(double));
#endif
#ifdef METALS
      SN_Z[z] = malloc(M_COUNT * sizeof(double));
#endif
#endif 
    }

  if(ThisTask == 0)
    {
      for (int z = 0; z < Z_COUNT; z++)
        {
          char zname[64];
          snprintf(zname, sizeof(zname), "Z=%g", Z_VALUES[z]);

          hid_t zgrp = my_H5Gopen(file_id, zname);

          for(int m = 0; m < M_COUNT; m++)
            {
              char mname[64];
              snprintf(mname, sizeof(mname), "M=%03d", (int)round((M_VALUES[m])));

              if (H5Lexists(zgrp, mname, H5P_DEFAULT) <= 0)
                {
                  terminate("Error loading stellar tables!");
                }

              hid_t mgrp = my_H5Gopen(zgrp, mname);
              
              hid_t d_age = my_H5Dopen(mgrp, "Age");
              hid_t d_frage = my_H5Dopen(mgrp, "FractionalAge");
              hid_t d_rad = my_H5Dopen(mgrp, "logRadius");
              hid_t d_tem = my_H5Dopen(mgrp, "logTemperature");

              hsize_t dims[1];
          
              hid_t space = H5Dget_space(d_age);
              H5Sget_simple_extent_dims(space, dims, NULL);
              H5Sclose(space);
              
              N[z][m] = (int)dims[0];

              Age[z][m] = malloc(N[z][m] * sizeof(double));
              FractionalAge[z][m] = malloc(N[z][m] * sizeof(double));
              logRadius[z][m] = malloc(N[z][m] * sizeof(double));
              logTemperature[z][m] = malloc(N[z][m] * sizeof(double));

              my_H5Dread(d_age, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, Age[z][m], "Age");          
              my_H5Dread(d_frage, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, FractionalAge[z][m], "FractionalAge");
              my_H5Dread(d_rad, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, logRadius[z][m], "logRadius");
              my_H5Dread(d_tem, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, logTemperature[z][m], "logTemperature");
              
              my_H5Dclose(d_age, "Age");
              my_H5Dclose(d_frage, "FractionalAge");
              my_H5Dclose(d_rad, "logRadius");
              my_H5Dclose(d_tem, "logTemperature");

#ifdef WINDS
              hid_t d_ml = my_H5Dopen(mgrp, "logMassLossRate");
#if GRACKLE_CHEMISTRY >= 1
              hid_t d_X  = my_H5Dopen(mgrp, "X");
              hid_t d_Y = my_H5Dopen(mgrp, "Y");
#endif
#ifdef METALS
              hid_t d_Z = my_H5Dopen(mgrp, "Z");
#endif
              hid_t d_wv = my_H5Dopen(mgrp, "logWindVelocity");

              logMassLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#if GRACKLE_CHEMISTRY >= 1
              WindX[z][m] = malloc(N[z][m] * sizeof(double));
              WindY[z][m] = malloc(N[z][m] * sizeof(double));
#endif
#ifdef METALS
              WindZ[z][m] = malloc(N[z][m] * sizeof(double));
#endif
              logWindVelocity[z][m] = malloc(N[z][m] * sizeof(double));

              my_H5Dread(d_ml, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, logMassLossRate[z][m], "logMassLossRate");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dread(d_X, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, WindX[z][m], "X");
              my_H5Dread(d_Y, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, WindY[z][m], "Y");
#endif
#ifdef METALS
              my_H5Dread(d_Z, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, WindZ[z][m], "Z");
#endif
              my_H5Dread(d_wv, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, logWindVelocity[z][m], "logWindVelocity");

              my_H5Dclose(d_ml, "logMassLossRate");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dclose(d_X, "X");
              my_H5Dclose(d_Y, "Y");
#endif
#ifdef METALS
              my_H5Dclose(d_Z, "Z");
#endif
              my_H5Dclose(d_wv, "logWindVelocity");
#endif

#ifdef STAR_RADIATION_ACTIVE
              hid_t d_energy = my_H5Dopen(mgrp, "logEnergy");
              hid_t d_photons = my_H5Dopen(mgrp, "logPhotons");

              double (*energy_buf)[WAVEBANDS] = malloc(N[z][m] * sizeof(*energy_buf));
              double (*photon_buf)[WAVEBANDS] = malloc(N[z][m] * sizeof(*photon_buf));

              my_H5Dread(d_energy, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, energy_buf, "logEnergy");
              my_H5Dread(d_photons, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, photon_buf, "logPhotons");

              my_H5Dclose(d_energy, "logEnergy");
              my_H5Dclose(d_photons, "logPhotons");

              for(int w = 0; w < WAVEBANDS; w++)
                {
                  logFlux[w][z][m] = malloc(N[z][m] * sizeof(WavebandData));
                  for(int i = 0; i < N[z][m]; i++)
                    {
                      logFlux[w][z][m][i].Energy = energy_buf[i][w]; 
                      logFlux[w][z][m][i].Photons = photon_buf[i][w]; 
                    }
                }

              free(energy_buf);
              free(photon_buf);
#endif

#ifdef SUPERNOVAE
              hid_t d_snml = my_H5Dopen(mgrp, "SN_MassLoss");
#if GRACKLE_CHEMISTRY >= 1
              hid_t d_snX = my_H5Dopen(mgrp, "SN_X");
              hid_t d_snY = my_H5Dopen(mgrp, "SN_Y");
#endif
#ifdef METALS
              hid_t d_snZ = my_H5Dopen(mgrp, "SN_Z");
#endif
              my_H5Dread(d_snml, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_MassLoss[z][m], "SN_MassLoss");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dread(d_snX, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_X[z][m], "SN_X");
              my_H5Dread(d_snY, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_Y[z][m], "SN_Y");
#endif
#ifdef METALS
              my_H5Dread(d_snZ, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL, H5P_DEFAULT, &SN_Z[z][m], "SN_Z");
#endif
              my_H5Dclose(d_snml, "SN_MassLoss");
#if GRACKLE_CHEMISTRY >= 1
              my_H5Dclose(d_snX, "SN_X");
              my_H5Dclose(d_snY, "SN_Y");
#endif
#ifdef METALS
              my_H5Dclose(d_snZ, "SN_Z");
#endif
#endif
              my_H5Gclose(mgrp, mname);
            }
          my_H5Gclose(zgrp, zname);
        }
      my_H5Fclose(file_id, filename);
    }

  for(int z = 0; z < Z_COUNT; z++)
    MPI_Bcast(N[z], M_COUNT, MPI_INT, 0, MPI_COMM_WORLD);

  for(int z = 0; z < Z_COUNT; z++)
    {
#ifdef SUPERNOVAE
      MPI_Bcast(SN_MassLoss[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if GRACKLE_CHEMISTRY >= 1
      MPI_Bcast(SN_X[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(SN_Y[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#ifdef METALS
      MPI_Bcast(SN_Z[z], M_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#endif
      for(int m = 0; m < M_COUNT; m++)
        if(N[z][m] > 0)
          {
            if(ThisTask != 0)
              {
                Age[z][m] = malloc(N[z][m] * sizeof(double));
                FractionalAge[z][m] = malloc(N[z][m] * sizeof(double));
                logRadius[z][m] = malloc(N[z][m] * sizeof(double));
                logTemperature[z][m] = malloc(N[z][m] * sizeof(double));

#ifdef WINDS
                logMassLossRate[z][m] = malloc(N[z][m] * sizeof(double));
#if GRACKLE_CHEMISTRY >= 1
                WindX[z][m] = malloc(N[z][m] * sizeof(double));
                WindY[z][m] = malloc(N[z][m] * sizeof(double));
#endif
#ifdef METALS
                WindZ[z][m] = malloc(N[z][m] * sizeof(double));
#endif
                logWindVelocity[z][m] = malloc(N[z][m] * sizeof(double));
#endif

#ifdef STAR_RADIATION_ACTIVE
                for(int w = 0; w < WAVEBANDS; w++)
                  logFlux[w][z][m] = malloc(N[z][m] * sizeof(WavebandData));
#endif
              }

            MPI_Bcast(Age[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(FractionalAge[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(logRadius[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(logTemperature[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);         

#ifdef WINDS
            MPI_Bcast(logMassLossRate[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if GRACKLE_CHEMISTRY >= 1
            MPI_Bcast(WindX[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(WindY[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
#ifdef METALS
            MPI_Bcast(WindZ[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
            MPI_Bcast(logWindVelocity[z][m], N[z][m], MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

#ifdef STAR_RADIATION_ACTIVE
            for(int w = 0; w < WAVEBANDS; w++)
              MPI_Bcast(logFlux[w][z][m], N[z][m] * sizeof(WavebandData), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
          }
    }
}