#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

static int bh_accretion_evaluate(int target, int mode, int threadid);
static void update_bh_accretion_rate(void);

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyDouble Vel[3];
  MyDouble AngularMomentum[3];
  MyFloat Hsml;
  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
  for(int j = 0; j < 3; j++)
    {
      in->Pos[j] = PPB(i).Pos[j];
      in->Vel[j] = PPB(i).Vel[j];
      in->AngularMomentum[j] = BhP[i].AngularMomentum[j];
    }
  in->Hsml = BhP[i].Hsml;
  in->Firstnode = firstnode;
}  

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
#ifdef BONDI_ACCRETION
  MyDouble VelocityGas[3];
  MyDouble VelocityGasCircular[3];
  MyDouble Density;
  MyDouble InternalEnergyGas;
#endif

#ifdef TORQUE_ACCRETION
  MyDouble TorqueMgas;
  MyDouble TorqueMstar;
  MyDouble TorqueMgasDisk;      /* Disk component gas mass */
  MyDouble TorqueMstarDisk;     /* Disk component stellar mass */ 
  MyDouble TorqueR0;
  MyDouble VelocityGasCircular[3];
#endif

#ifdef ADP_ACCRETION
  MyDouble ADP_CapturedMass;
  MyDouble ADP_ReservoirMass;
  MyDouble ADP_DiscMass;
  MyDouble VelocityGas[3];
  MyDouble VelocityGasCircular[3];
#endif

//#ifdef INFALL_ACCRETION
//  MyDouble Accretion;
//#endif
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *  particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *  communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
#ifdef BONDI_ACCRETION
      for(int j = 0; j < 3; j++)
        {
          BhP[i].VelocityGas[j] = out->VelocityGas[j];
          BhP[i].VelocityGasCircular[j] = out->VelocityGasCircular[j];
        }
      BhP[i].Density = out->Density;
      BhP[i].InternalEnergyGas = out->InternalEnergyGas;
#endif

#ifdef TORQUE_ACCRETION
      BhP[i].TorqueMgas = out->TorqueMgas;
      BhP[i].TorqueMstar = out->TorqueMstar;
      BhP[i].TorqueMgasDisk = out->TorqueMgasDisk;
      BhP[i].TorqueMstarDisk = out->TorqueMstarDisk;
      BhP[i].TorqueR0 = out->TorqueR0;
      for(int j = 0; j < 3; j++)
        BhP[i].VelocityGasCircular[j] = out->VelocityGasCircular[j];
#endif

#ifdef ADP_ACCRETION
      BhP[i].ADP_CapturedMass = out->ADP_CapturedMass;
      BhP[i].ADP_ReservoirMass = out->ADP_ReservoirMass;
      BhP[i].ADP_DiscMass = out->ADP_DiscMass;
      for(int j = 0; j < 3; j++)
        {
          BhP[i].VelocityGas[j] = out->VelocityGas[j];
          BhP[i].VelocityGasCircular[j] = out->VelocityGasCircular[j];
        }
#endif

//#ifdef INFALL_ACCRETION
//      BhP[i].Accretion = out->Accretion;
//#endif
    }
  else /* combine */
    {
#ifdef BONDI_ACCRETION
      for(int j = 0; j < 3; j++)
        {
          BhP[i].VelocityGas[j] += out->VelocityGas[j];
          BhP[i].VelocityGasCircular[j] += out->VelocityGasCircular[j];
        }
      BhP[i].Density += out->Density;
      BhP[i].InternalEnergyGas += out->InternalEnergyGas;
#endif

#ifdef TORQUE_ACCRETION
      BhP[i].TorqueMgas += out->TorqueMgas;
      BhP[i].TorqueMstar += out->TorqueMstar;
      BhP[i].TorqueMgasDisk += out->TorqueMgasDisk;
      BhP[i].TorqueMstarDisk += out->TorqueMstarDisk;
      if(out->TorqueR0 > BhP[i].TorqueR0)
        BhP[i].TorqueR0 = out->TorqueR0;
      for(int j = 0; j < 3; j++)
        BhP[i].VelocityGasCircular[j] += out->VelocityGasCircular[j];
#endif

#ifdef ADP_ACCRETION
      BhP[i].ADP_CapturedMass += out->ADP_CapturedMass;
      BhP[i].ADP_ReservoirMass += out->ADP_ReservoirMass;
      BhP[i].ADP_DiscMass += out->ADP_DiscMass;
      for(int j = 0; j < 3; j++)
        {
          BhP[i].VelocityGas[j] += out->VelocityGas[j];
          BhP[i].VelocityGasCircular[j] += out->VelocityGasCircular[j];
        }
#endif

//#ifdef INFALL_ACCRETION
//      BhP[i].Accretion += out->Accretion; 
//#endif
    }
}


#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i, idx;
  int j, threadid = get_thread_num();

  for(j = 0; j < NTask; j++)
    Thread[threadid].Exportflag[j] = -1;

  while(1)
    {
      if(Thread[threadid].ExportSpace < MinSpace)
        break;

      //i = NextParticle++;

      //if(i >= NumBhs)
      //break;
        
      idx = NextParticle++;

      if(idx >= TimeBinsBh.NActiveParticles)
        break;

      i = TimeBinsBh.ActiveParticleList[idx];
        
      bh_accretion_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
    }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  int threadid = get_thread_num();

  while(1)
    {
      i = cnt++;

      if(i >= Nimport)
        break;

      bh_accretion_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
    }
}

/*! \brief Main function of SPH density calculation.
 *
 *  This function computes the local density for each active SPH particle and
 *  the number of weighted neighbors in the current smoothing radius. If a
 *  particle with its smoothing region is fully inside the local domain, it is
 *  not exported to the other processors. The function also detects particles
 *  that have a number of neighbors outside the allowed tolerance range. For
 *  these particles, the smoothing length is adjusted accordingly, and the
 *  computation is called again.
 *
 *  \return void
 */
void bh_accretion(void)
{
  int i, idx;

  CPU_Step[CPU_MISC] += measure_time();

  generic_set_MaxNexport();

  generic_comm_pattern(TimeBinsBh.NActiveParticles, kernel_local, kernel_imported);

#ifdef TORQUE_ACCRETION
  for(idx = 0; idx < TimeBinsBh.NActiveParticles; idx++)
    {
      i = TimeBinsBh.ActiveParticleList[idx];
      MyDouble M_total = BhP[i].TorqueMgas; //+ BhP[i].TorqueMstar; TODO
      MyDouble M_disk = BhP[i].TorqueMgasDisk; //+ BhP[i].TorqueMstarDisk; 
      BhP[i].TorqueFd = (M_total > 0) ? M_disk / M_total : 0.0;
    }
#endif

  update_bh_accretion_rate();

  CPU_Step[CPU_INIT] += measure_time();
}

/*! \brief Inner function of the SPH density calculation
 *
 *  This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 *
 *  \param[in] target Index of particle in local data/import buffer.
 *  \param[in] mode Mode in which function is called (local or impored data).
 *  \param[in] threadid ID of local thread.
 *
 *  \return 0
 */
static int bh_accretion_evaluate(int target, int mode, int threadid)
{
  int i, n, numnodes, *firstnode; 
  double h, h2, r, r2, wk;
  double dx, dy, dz, dvx, dvy, dvz; 
  MyDouble *pos, *vel, *angular_momentum;
  
  data_in local, *target_data;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  vel = target_data->Vel;
  angular_momentum = target_data->AngularMomentum;
  h = target_data->Hsml;

#ifdef BONDI_ACCRETION
  MyDouble velocity_gas[3], velocity_gas_circular[3];
  for(int j = 0; j < 3; j++)
    velocity_gas[j] = velocity_gas_circular[j] = 0;
  MyDouble density, internal_energy_gas;
  density = internal_energy_gas = 0;
#endif

#ifdef TORQUE_ACCRETION
  /* Torque accretion model from Angles-Alcazar et al. 2016 */
  MyDouble torque_Mgas = 0.0;  /* Total gas mass within R0 */
  MyDouble torque_Mstar = 0.0;  /* Total stellar mass within R0 */
  MyDouble torque_Mgas_disk = 0.0;  /* Disk component gas mass */
  MyDouble torque_Mstar_disk = 0.0;  /* Disk component stellar mass */
  MyDouble R0_torque = All.TorqueR0;   // 0.2-0.3 kpc in code units
  MyDouble R0_torque2 = R0_torque * R0_torque;
  /* For angular momentum accretion tracking */
  MyDouble velocity_gas_circular[3];
  for(int j = 0; j < 3; j++)
    velocity_gas_circular[j] = 0;
#endif

#ifdef ADP_ACCRETION
  MyDouble adp_captured_mass = 0.0;
  MyDouble Racc = h * 24 ;
  MyDouble Racc2 = Racc * Racc;
#endif

//#ifdef INFALL_ACCRETION
//  MyDouble accretion = 0;
//  double rbh = h;
//  double rbh2 = rbh * rbh;
//#endif

  double hinv, hinv3, hinv4, u, dwk;

  h2 = h * h;
  hinv = 1.0 / h;
#ifndef TWODIMS
  hinv3 = hinv * hinv * hinv;
#else  /* #ifndef  TWODIMS */
  hinv3 = hinv * hinv / boxSize_Z;
#endif /* #ifndef  TWODIMS #else */
  hinv4 = hinv3 * hinv;

  int nfound = ngb_treefind_variable_threads(pos, h, target, mode, threadid, numnodes, firstnode);

  for(n = 0; n < nfound; n++)
    {
      i = Thread[threadid].Ngblist[n];

/* compute bh->cell position vectors: posBhP-posSphP */
      dx = pos[0] - P[i].Pos[0];
      dy = pos[1] - P[i].Pos[1];
      dz = pos[2] - P[i].Pos[2];

/* compute bh->cell velocity vectors: posBhP-posSphP */
      dvx = vel[0] - P[i].Vel[0];
      dvy = vel[1] - P[i].Vel[1];
      dvz = vel[2] - P[i].Vel[2];

#ifndef REFLECTIVE_X
      if(dx > boxHalf_X)
        dx -= boxSize_X;
      if(dx < -boxHalf_X)
        dx += boxSize_X;
#endif /* #ifndef REFLECTIVE_X */

#ifndef REFLECTIVE_Y
      if(dy > boxHalf_Y)
        dy -= boxSize_Y;
      if(dy < -boxHalf_Y)
        dy += boxSize_Y;
#endif /* #ifndef REFLECTIVE_Y */

#ifndef REFLECTIVE_Z
      if(dz > boxHalf_Z)
        dz -= boxSize_Z;
      if(dz < -boxHalf_Z)
        dz += boxSize_Z;
#endif /* #ifndef REFLECTIVE_Z */
      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 < h2)
        {
          r = sqrt(r2);

          u = r * hinv;

          kernel(u, hinv3, hinv4, &wk, &dwk);

          double v_cross[3];
          v_cross[0] = dy * dvz - dz * dvy;
          v_cross[1] = dz * dvx - dx * dvz;
          v_cross[2] = dx * dvy - dy * dvx;

          double rho = (SphP[i].Density > 0) ? SphP[i].Density : 1;

#ifdef BONDI_ACCRETION
          /* comute relative velocities, 
          relative specific angular momenta and internal energy of gas */

          velocity_gas[0] += dvx * P[i].Mass / rho * wk;
          velocity_gas[1] += dvy * P[i].Mass / rho * wk;
          velocity_gas[2] += dvz * P[i].Mass / rho * wk;

          velocity_gas_circular[0] += v_cross[0] * P[i].Mass / rho * wk;
          velocity_gas_circular[1] += v_cross[1] * P[i].Mass / rho * wk;
          velocity_gas_circular[2] += v_cross[2] * P[i].Mass / rho * wk;

          density += P[i].Mass * wk;

          internal_energy_gas += SphP[i].Utherm * P[i].Mass / rho * wk;
#endif

#ifdef TORQUE_ACCRETION
          if(r2 < R0_torque2)
            {
              double v_phi = (v_cross[0]*angular_momentum[0] + v_cross[1]*angular_momentum[1] + v_cross[2]*angular_momentum[2]) / r;
          
              int is_disk = (v_phi > 0.0);
          
              /* Accumulate total masses */
              if(P[i].Type == 0)  /* Gas */
                {
                  torque_Mgas += P[i].Mass; 
                  
                  if(is_disk)
                    torque_Mgas_disk += P[i].Mass; 
              
                  /* Accumulate circular velocity for angular momentum tracking */
                  velocity_gas_circular[0] += v_cross[0] * P[i].Mass / rho * wk;
                  velocity_gas_circular[1] += v_cross[1] * P[i].Mass / rho * wk;
                  velocity_gas_circular[2] += v_cross[2] * P[i].Mass / rho * wk;
                }
          
            /*else if(P[i].Type == 4) TODO
            {
              torque_Mstar += P[i].Mass; 
              if(is_disk)
                torque_Mstar_disk += P[i].Mass; 
            }*/
            
            }
#endif 

#ifdef ADP_ACCRETION
          if(r2 < Racc2)
            adp_captured_mass += P[i].Mass;
#endif

//#ifdef INFALL_ACCRETION
          /* cell nibbled */
//          if(r < 2*rbh) 
//            {
//              accretion += P[j].Mass * exp(-r2/(2*rbh2));
//              P[j].Mass -= P[j].Mass * exp(-r2/(2*rbh2));  
//            }
//#endif

        } // if(r2 < h2)
    } // for(n = 0; n < nfound; n++)

#ifdef BONDI_ACCRETION
  for(int j = 0; j < 3; j++)
    {
      out.VelocityGas[j] = velocity_gas[j];
      out.VelocityGasCircular[j] = velocity_gas_circular[j];
    }
  out.Density = density;
  out.InternalEnergyGas = internal_energy_gas;
#endif

#ifdef TORQUE_ACCRETION
  out.TorqueMgas = torque_Mgas;
  out.TorqueMstar = torque_Mstar;
  out.TorqueMgasDisk = torque_Mgas_disk;
  out.TorqueMstarDisk = torque_Mstar_disk;
  out.TorqueR0 = h; 
  for(int j = 0; j < 3; j++) 
    out.VelocityGasCircular[j] = velocity_gas_circular[j];
#endif

#ifdef ADP_ACCRETION
  out.ADP_CapturedMass = adp_captured_mass;
#endif

//#ifdef INFALL_ACCRETION
//  out.Accretion = accretion;
//#endif

  /* now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

#ifdef BONDI_ACCRETION 
static void update_bh_accretion_rate(void)
{
  /* calculate bondi accretion rate */
  int i;
  double density, pressure, sound_speed, velocity_gas_norm;
  double denominator, denominator_inv, BondiRate, EddingtonRate;
  double accretion_rate, acc_max, acc_rate_for_print;

  accretion_rate = acc_max = acc_rate_for_print = 0;

  for(i = 0; i < NumBhs; i++)
    {
      /* get pressure */
      if(BhP[i].Density > 0)
        {  
          density = BhP[i].Density;
          pressure = GAMMA_MINUS1 * density * BhP[i].InternalEnergyGas;

          /* get soundspeed */
          sound_speed = sqrt(GAMMA * pressure / density);
      
          velocity_gas_norm = sqrt(BhP[i].VelocityGas[0]*BhP[i].VelocityGas[0] + 
          BhP[i].VelocityGas[1]*BhP[i].VelocityGas[1] + BhP[i].VelocityGas[2]*BhP[i].VelocityGas[2]); 

          denominator = (sound_speed*sound_speed + velocity_gas_norm*velocity_gas_norm);
          if(denominator > 0)
            {
              denominator_inv = 1. / sqrt(denominator);
              BondiRate = 4. * M_PI * All.G * All.G * PPB(i).Mass * PPB(i).Mass * density *
              denominator_inv * denominator_inv * denominator_inv;
            }
          else
            terminate("Invalid denominator in Bondi Accretion Rate");
        }
      else
        BondiRate = 0;
  
      /* limit by Eddington accretion rate */
      EddingtonRate = 4. * M_PI * GRAVITY * (PPB(i).Mass * All.UnitMass_in_g) * PROTONMASS / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *=  (All.UnitTime_in_s / All.UnitMass_in_g);
      accretion_rate = fmin(BondiRate, EddingtonRate);
      
      BhP[i].AccretionRate  = accretion_rate;

      /* Track maximum for output */
      if(accretion_rate > acc_max)
        acc_max = accretion_rate;
    }
 
  MPI_Allreduce(&acc_max, &acc_rate_for_print, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  mpi_printf("BLACK_HOLES: Black hole max Bondi accretion rate: %e (code units)\n", acc_rate_for_print);
}
#endif

#ifdef TORQUE_ACCRETION
static void update_bh_accretion_rate(void)
{
  /* Calculate Torque-limited accretion rate */
  int i;
  double M_BH, M_gas, M_star, M_enc, M_gas_disk, M_star_disk, M_disk;
  double R0, f_d, f_gas, f0;
  double torque_rate, EddingtonRate;
  double accretion_rate, acc_max, acc_rate_for_print;

  accretion_rate = acc_max = acc_rate_for_print = 0;
  
  for(i = 0; i < NumBhs; i++)
    {
      M_BH        = PPB(i).Mass;
      M_gas       = BhP[i].TorqueMgas;       // Total gas mass within R0
      M_star      = BhP[i].TorqueMstar;      // Total stellar mass within R0
      M_gas_disk  = BhP[i].TorqueMgasDisk;   // Disk component gas mass 
      M_star_disk = BhP[i].TorqueMstarDisk;  // Disk component stellar mass 
      R0          = BhP[i].TorqueR0;         // Aperture radius
      f_d         = BhP[i].TorqueFd;         // Disk fraction

      if(R0 <= 0 || (M_gas + M_star) <= 0)
        {
          BhP[i].AccretionRate = 0;
          continue;
        }

      /* Disk mass (gas + stars in disk component) */
      M_disk = M_gas_disk + M_star_disk;
      
      /* Total enclosed mass */
      M_enc = M_gas + M_star;
      
      /* Recompute disk fraction  */
      if(M_enc > 0)
        f_d = M_disk / M_enc;
      else
        f_d = 0.0;
        
      /* Gas fraction within disk */
      if(M_disk > 0)
        f_gas = M_gas_disk / M_disk;  // Only disk gas / disk mass
      else
        f_gas = 1.0;  // Default to pure gas if no disk

      /* Ensure physical values */
      if(f_d < 0.0) f_d = 0.0;
      if(f_d > 1.0) f_d = 1.0;
      if(f_gas < 0.0) f_gas = 0.0;
      if(f_gas > 1.0) f_gas = 1.0;
      
      /* If no disk, no torque-driven accretion */
      if(f_d < 1e-6 || M_disk < 1e-6)
        {
          BhP[i].AccretionRate = 0;
          continue;
        }

      /* Computing f0 */
      /* f0 ≈ 0.31 * f_d^2 * (M_disk / 10^9 Msun)^(-1/3) */
      f0 = 0.31 * f_d * f_d * pow(M_disk * All.UnitMass_in_g / 1e9 / SOLAR_MASS, -1.0/3.0);

      /* Suppression factor */
      double suppression = 1.0;
      if(f_gas > 0)
        suppression = 1.0 / (1.0 + f0 / f_gas);
      else
        suppression = 0.0;  // No gas, no accretion

      /* Torque-limited accretion rate (Angles-Alcazar et al. 2016, Equation 2) */

      torque_rate = All.Epsilon_T                                      /* Normalization */
                  * pow(f_d, 2.5)                                      /* f_d^(5/2) */
                  * pow(M_BH * All.UnitMass_in_g / 1e8 / SOLAR_MASS, 1.0/6.0)  /* (M_BH/1e8 Msun)^(1/6) */
                  * (M_disk * All.UnitMass_in_g / 1e9 / SOLAR_MASS)   /* (M_disk/1e9 Msun) */
                  * pow(R0 * All.UnitLength_in_cm / (100.0 * PARSEC), -1.5)    /* (R0/100pc)^(-3/2) */
                  * suppression;                                       /* 1/(1 + f0/f_gas) */

      /* Convert from Msun/yr to code units (code mass / code time) */
      torque_rate *= SOLAR_MASS / All.UnitMass_in_g;         
      torque_rate *= All.UnitTime_in_s / SEC_PER_YEAR;       

      /* Eddington limit (can allow up to 10× super-Eddington as in paper) */
      EddingtonRate = 4.0 * M_PI * GRAVITY * (M_BH * All.UnitMass_in_g) * PROTONMASS
                      / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *= (All.UnitTime_in_s / All.UnitMass_in_g);
      
      /* Allow up to 10× Eddington (Angles-Alcazar et al. 2016, Section 2.3) */
      double max_accretion = 10.0 * EddingtonRate;

      /* Apply Eddington limit */
      accretion_rate = fmin(torque_rate, max_accretion);

      /* Store the accretion rate */
      BhP[i].AccretionRate = accretion_rate;
      
      /* Track maximum for output */
      if(accretion_rate > acc_max)
        acc_max = accretion_rate;
    }
 
  MPI_Allreduce(&acc_max, &acc_rate_for_print, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  mpi_printf("BLACK_HOLES: Black hole max torque-limited accretion rate: %e (code units)\n", acc_rate_for_print);
}
#endif

#ifdef ADP_ACCRETION
static void update_bh_accretion_rate(void)
{
  int i;
  int bin;
  double dt;
  double M_BH;
  double Mcap;
  double M_res;
  double M_disc;
  double dM_to_disc, mdot_visc, mdot_cap, dM_bh;
  double EddingtonRate; 
  double accretion_rate, acc_max, acc_rate_for_print;

  accretion_rate = acc_max = acc_rate_for_print = 0;

  for(i = 0; i < NumBhs; i++)
    {
      bin = BhP[i].TimeBinBh;
      dt = (bin ? (((integertime)1) << bin) : 0) * All.Timebase_interval;

      if(dt <= 0) continue;

      M_BH = PPB(i).Mass;

      Mcap = BhP[i].ADP_CapturedMass;   /* set by bh_density this step    */
      M_res = BhP[i].ADP_ReservoirMass;  /* carried over from previous step */
      M_disc = BhP[i].ADP_DiscMass;       /* carried over from previous step */

      if(Mcap < 0) Mcap = 0;

      /* 
       * Stage 1 → 2: Captured mass enters the reservoir immediately.
           Gas that crossed Racc is not yet on the disc — it still has angular
           momentum and must circularise first.  The reservoir drains into the
           disc on the dynamical / capture timescale ADP_tcap. */
      M_res += Mcap;
      BhP[i].ADP_CapturedMass = 0;

      /* How much flows from reservoir → disc this timestep?
         dM = M_res * (dt / tcap).
         If ADP_tcap == 0 (instantaneous) dump everything at once. */
      if(All.ADP_tcap > 0)
        dM_to_disc = M_res * (dt / All.ADP_tcap);
      else
        dM_to_disc = M_res;   /* instantaneous: reservoir empties each step */

      if(dM_to_disc > M_res) dM_to_disc = M_res;
      if(dM_to_disc < 0)     dM_to_disc = 0;

      M_res  -= dM_to_disc;
      M_disc += dM_to_disc;

      /*  Stage 2 → 3: Disc drains onto the BH on the viscous timescale.
          Mdot_BH = min( Mdisc / tvisc , Mdot_Edd ) */
      if(All.ADP_tvisc > 0)
        mdot_visc = M_disc / All.ADP_tvisc;
      else
        mdot_visc = (dt > 0) ? M_disc / dt : 0;   /* fallback: drain in one step */

      if(mdot_visc < 0) mdot_visc = 0;

      EddingtonRate = 4.0 * M_PI * GRAVITY * (M_BH * All.UnitMass_in_g) * PROTONMASS
                    / (All.Epsilon_r * CLIGHT * THOMPSON);
      EddingtonRate *= (All.UnitTime_in_s / All.UnitMass_in_g);

      mdot_cap = All.ADP_EddFactor * EddingtonRate;

      accretion_rate = fmin(mdot_visc, mdot_cap);

      dM_bh = accretion_rate * dt;

      if(dM_bh > M_disc)
        {
          dM_bh          = M_disc;
          accretion_rate = (dt > 0) ? dM_bh / dt : 0;
        }
      if(dM_bh < 0) dM_bh = 0;

      M_disc -= dM_bh;

      BhP[i].ADP_ReservoirMass = M_res;
      BhP[i].ADP_DiscMass = M_disc;
      BhP[i].AccretionRate = accretion_rate;

      //debug
      //mpi_printf("ADP BH %d: Mcap=%e  Mres=%e  Mdisc=%e  mdot_visc=%e  mdot_Edd=%e  Mdot_BH=%e\n",
      //           i, Mcap, M_res, M_disc, mdot_visc, EddingtonRate, accretion_rate);

      /* Track maximum for output */
      if(accretion_rate > acc_max)
        acc_max = accretion_rate;
    }
 
  MPI_Allreduce(&acc_max, &acc_rate_for_print, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); 
  mpi_printf("BLACK_HOLES: Black hole max ADP accretion rate: %e (code units)\n", acc_rate_for_print);
}
#endif