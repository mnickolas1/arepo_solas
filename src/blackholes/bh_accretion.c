#include <stdlib.h>       
#include <math.h>
#include <gsl/gsl_math.h>              
#include <mpi.h>            
  
#include "../main/allvars.h"
#include "../main/proto.h"

#include "../domain/domain.h"

static int bh_accretion_evaluate(int target, int mode, int threadid);

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
  int j, n, numnodes, *firstnode; 
  int numngb; 
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
  MyDouble internal_energy_gas = 0;
  MyDouble velocity_gas[3], velocity_gas_circular[3];
  for(int j = 0; j < 3; j++)
    velocity_gas[j] = velocity_gas_circular[j] = 0;
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
  MyDouble velocity_gas[3], velocity_gas_circular[3];
  for(int j = 0; j < 3; j++)
    velocity_gas[j] = velocity_gas_circular[j] = 0;
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
          numngb++;

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

          velocity_gas[0] += dvx * P[i].Mass / rho*wk;
          velocity_gas[1] += dvy * P[i].Mass / rho*wk;
          velocity_gas[2] += dvz * P[i].Mass / rho*wk;

          velocity_gas_circular[0] -= v_cross[0] * P[i].Mass / rho*wk;
          velocity_gas_circular[1] -= v_cross[1] * P[i].Mass / rho*wk;
          velocity_gas_circular[2] -= v_cross[2] * P[i].Mass / rho*wk;

          internal_energy_gas += SphP[i].Utherm * P[i].Mass / rho*wk;
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
                  
                  velocity_gas_circular[0] += v_cross[0] * P[i].Mass / rho;
                  velocity_gas_circular[1] += v_cross[1] * P[i].Mass / rho;
                  velocity_gas_circular[2] += v_cross[2] * P[i].Mass / rho;
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
          
          double v_phi = (v_cross[0]*angular_momentum[0] + v_cross[1]*angular_momentum[1] + v_cross[2]*angular_momentum[2]) / r;
          
          int is_disk = (v_phi > 0.0);
          
          /* Accumulate total masses */
          if(P[j].Type == 0)  /* Gas */
            {
              adp_captured_mass += mass_j; 
              if(is_disk)
                adp_captured_mass += mass_j; 
              
              /* Accumulate circular velocity for angular momentum tracking */
                  
              velocity_gas_circular[0] += v_cross[0] * P[i].Mass / rho;
              velocity_gas_circular[1] += v_cross[1] * P[i].Mass / rho;
              velocity_gas_circular[2] += v_cross[2] * P[i].Mass / rho;
            }
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
  for(int j = 0; j < 3; j++)
    {
      out.VelocityGas[j] = velocity_gas[j];
      out.VelocityGasCircular[j] = velocity_gas_circular[j];
    }
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