/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/domain_exchange.c
 * \date        05/2018
 * \brief       Algorithms for exchanging particle data and associated
 *              rearrangements.
 * \details     This includes changing the size of the P and SphP arrays as
 *              well as the particle exchange routine itself.
 *              contains functions:
 *                void domain_resize_storage(int count_get, int count_get_sph,
 *                  int option_flag)
 *                void domain_exchange(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 05.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"
#include "domain.h"

/*! \brief Changes memory allocation if necessary for particle and cell data.
 *
 *  If the memory usage due to a net import or export of particles changes
 *  above a certain tolerance, the P and SphP structures need to be
 *  reallocated.
 *
 *  \param[in] count get How many particles are imported?
 *  \param[in] count_get_sph How many cells are imported?
 *  \param[in] option_flag Options for reallocating peanokey or ngbtree.
 *
 *  \return void
 */
void domain_resize_storage(int count_get, int count_get_sph, int option_flag)
{
  int load        = NumPart + count_get;
  int sphload     = NumGas + count_get_sph;
  int loc_data[2] = {load, sphload}, res[2];

  MPI_Allreduce(loc_data, res, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int max_load    = res[0];
  int max_sphload = res[1];

  if(max_load > (1.0 - ALLOC_TOLERANCE) * All.MaxPart || max_load < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPart)
    {
      All.MaxPart = max_load / (1.0 - 2 * ALLOC_TOLERANCE);
      reallocate_memory_maxpart();

      if(option_flag == 1)
        Key = (peanokey *)myrealloc_movable(Key, sizeof(peanokey) * All.MaxPart);
    }

  if(max_sphload >= (1.0 - ALLOC_TOLERANCE) * All.MaxPartSph || max_sphload < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPartSph)
    {
      All.MaxPartSph = max_sphload / (1.0 - 2 * ALLOC_TOLERANCE);
      if(option_flag == 2)
        {
          if(All.MaxPartSph > Ngb_MaxPart)
            ngb_treemodifylength(All.MaxPartSph - Ngb_MaxPart);
        }
      reallocate_memory_maxpartsph();
    }
}

#ifdef BLACKHOLES
void domain_resize_storage_bhs(int count_get_bhs)
{
  /*int bhload       = NumBhs + count_get_bhs;
  int loc_data_bh  = bhload;
  int res_bh;

  MPI_Allreduce(&loc_data_bh, &res_bh, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int max_bhload   = res_bh;
  
  if(max_bhload > (1.0 - ALLOC_TOLERANCE) * All.MaxPartBhs || max_bhload < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPartBhs)
    {
      All.MaxPartBhs = max_bhload / (1.0 - 2 * ALLOC_TOLERANCE);
      reallocate_memory_maxpartbhs();
    }*/
}
#endif

#ifdef STARS
void domain_resize_storage_stars(int count_get_stars)
{
  int starload       = NumStars + count_get_stars;
  int loc_data_star  = starload;
  int res_star;

  MPI_Allreduce(&loc_data_star, &res_star, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int max_starload   = res_star;

  if(max_starload > 250)
  {
  
  if(max_starload > (1.0 - ALLOC_TOLERANCE) * All.MaxPartStars || max_starload < (1.0 - 3 * ALLOC_TOLERANCE) * All.MaxPartStars)
    {
      All.MaxPartStars = max_starload / (1.0 - 2 * ALLOC_TOLERANCE);
      reallocate_memory_maxpartstars();
    }

  }
}
#endif

/*! \brief Exchanges particles and cells according to new domain decomposition.
 *
 *  Communicates particles and cells to their new task. P and SphP arrays are
 *  changed in size accordingly.
 *
 *  \return void
 */
void domain_exchange(void)
{
  double t0 = second();

  int count_togo = 0, count_togo_sph = 0, count_get = 0, count_get_sph = 0;
  int *count, *count_sph, *offset, *offset_sph;
  int *count_recv, *count_recv_sph, *offset_recv, *offset_recv_sph;
  int i, n, no, target;
  struct particle_data *partBuf;
  struct sph_particle_data *sphBuf;
#ifdef BLACKHOLES
  int count_togo_bhs = 0, count_get_bhs = 0;
  int *count_bhs, *offset_bhs;
  int *count_recv_bhs, *offset_recv_bhs;
  struct bh_particle_data *bhBuf;
#endif
#ifdef STARS
  int count_togo_stars = 0, count_get_stars = 0;
  int *count_stars, *offset_stars;
  int *count_recv_stars, *offset_recv_stars;
  struct star_particle_data *sBuf;
#endif

  peanokey *keyBuf;

  long long sumtogo = 0;

  for(i = 0; i < NTask; i++)
    sumtogo += toGo[i];

  sumup_longs(1, &sumtogo, &sumtogo);

  count           = (int *)mymalloc_movable(&count, "count", NTask * sizeof(int));
  count_sph       = (int *)mymalloc_movable(&count_sph, "count_sph", NTask * sizeof(int));
  offset          = (int *)mymalloc_movable(&offset, "offset", NTask * sizeof(int));
  offset_sph      = (int *)mymalloc_movable(&offset_sph, "offset_sph", NTask * sizeof(int));
  count_recv      = (int *)mymalloc_movable(&count_recv, "count_recv", NTask * sizeof(int));
  count_recv_sph  = (int *)mymalloc_movable(&count_recv_sph, "count_recv_sph", NTask * sizeof(int));
  offset_recv     = (int *)mymalloc_movable(&offset_recv, "offset_recv", NTask * sizeof(int));
  offset_recv_sph = (int *)mymalloc_movable(&offset_recv_sph, "offset_recv_sph", NTask * sizeof(int));
#ifdef BLACKHOLES
  count_bhs        = (int *)mymalloc_movable(&count_bhs, "count_bhs", NTask * sizeof(int));
  offset_bhs       = (int *)mymalloc_movable(&offset_bhs, "offset_bhs", NTask * sizeof(int));
  count_recv_bhs   = (int *)mymalloc_movable(&count_recv_bhs, "count_recv_bhs", NTask * sizeof(int));
  offset_recv_bhs  = (int *)mymalloc_movable(&offset_recv_bhs, "offset_recv_bhs", NTask * sizeof(int));
#endif
#ifdef STARS
  count_stars        = (int *)mymalloc_movable(&count_stars, "count_stars", NTask * sizeof(int));
  offset_stars       = (int *)mymalloc_movable(&offset_stars, "offset_stars", NTask * sizeof(int));
  count_recv_stars   = (int *)mymalloc_movable(&count_recv_stars, "count_recv_stars", NTask * sizeof(int));
  offset_recv_stars  = (int *)mymalloc_movable(&offset_recv_stars, "offset_recv_stars", NTask * sizeof(int));
#endif

  int prec_offset;
  int *decrease;

  decrease = (int *)mymalloc_movable(&decrease, "decrease", NTask * sizeof(int));

  for(i = 1, offset_sph[0] = 0, decrease[0] = 0; i < NTask; i++)
    {
      offset_sph[i] = offset_sph[i - 1] + toGoSph[i - 1];
      decrease[i]   = toGoSph[i - 1];
    }
#ifdef BLACKHOLES
  for(i = 1, offset_bhs[0] = 0; i < NTask; i++)
    {
      offset_bhs[i] = offset_bhs[i - 1] + toGoBhs[i - 1];
    }
#endif
#ifdef STARS
  for(i = 1, offset_stars[0] = 0; i < NTask; i++)
    {
      offset_stars[i] = offset_stars[i - 1] + toGoStars[i - 1];
    }
#endif

  prec_offset = offset_sph[NTask - 1] + toGoSph[NTask - 1];

  offset[0] = prec_offset;
  for(i = 1; i < NTask; i++)
    offset[i] = offset[i - 1] + (toGo[i - 1] - decrease[i]);

  myfree(decrease);

  for(i = 0; i < NTask; i++)
    {
      count_togo += toGo[i];
      count_togo_sph += toGoSph[i];
      count_get += toGet[i];
      count_get_sph += toGetSph[i];
#ifdef BLACKHOLES
      count_togo_bhs += toGoBhs[i];
      count_get_bhs += toGetBhs[i];
#endif
#ifdef STARS
      count_togo_stars += toGoStars[i];
      count_get_stars += toGetStars[i];
#endif
    }

  partBuf = (struct particle_data *)mymalloc_movable(&partBuf, "partBuf", count_togo * sizeof(struct particle_data));
  sphBuf  = (struct sph_particle_data *)mymalloc_movable(&sphBuf, "sphBuf", count_togo_sph * sizeof(struct sph_particle_data));
#ifdef BLACKHOLES
  bhBuf = (struct bh_particle_data *)mymalloc_movable(&bhBuf, "bhBuf", count_togo_bhs * sizeof(struct bh_particle_data));
#endif
#ifdef STARS
  sBuf = (struct star_particle_data *)mymalloc_movable(&sBuf, "sBuf", count_togo_stars * sizeof(struct star_particle_data));
#endif
  
  keyBuf = (peanokey *)mymalloc_movable(&keyBuf, "keyBuf", count_togo * sizeof(peanokey));

  for(i = 0; i < NTask; i++)
    {
      count[i] = count_sph[i] = 0;
#ifdef BLACKHOLES
      count_bhs[i] = 0;
#endif
#ifdef STARS
      count_stars[i] = 0;
#endif
    }

  for(n = 0; n < NumPart; n++)
    {
      no = 0;

      peanokey mask = ((peanokey)7) << (3 * (BITS_PER_DIMENSION - 1));
      int shift     = 3 * (BITS_PER_DIMENSION - 1);

      while(topNodes[no].Daughter >= 0)
        {
          no = topNodes[no].Daughter + (int)((Key[n] & mask) >> shift);
          mask >>= 3;
          shift -= 3;
        }

      no = topNodes[no].Leaf;

      target = DomainTask[no];

      if(target != ThisTask)
        {
          /* copy this particle into the exchange buffer */
          if(P[n].Type == 0)
            {
              partBuf[offset_sph[target] + count_sph[target]] = P[n];
              keyBuf[offset_sph[target] + count_sph[target]]  = Key[n];
              sphBuf[offset_sph[target] + count_sph[target]]  = SphP[n];
              count_sph[target]++;
            }
#ifdef BLACKHOLES
          else if(P[n].Type == 5)
            {
              bhBuf[offset_bhs[target] + count_bhs[target]] = BPP(n);
              partBuf[offset[target] + count[target]] = P[n];
              keyBuf[offset[target] + count[target]]  = Key[n];
              
              count_bhs[target]++;
              count[target]++;  
            }
#endif
#ifdef STARS
          else if(P[n].Type == 4)
            {
              sBuf[offset_stars[target] + count_stars[target]] = SPP(n);
              partBuf[offset[target] + count[target]] = P[n];
              keyBuf[offset[target] + count[target]]  = Key[n];
              
              count_stars[target]++;
              count[target]++;  
            }
#endif
          else
            {
              partBuf[offset[target] + count[target]] = P[n];
              keyBuf[offset[target] + count[target]]  = Key[n];
              count[target]++;
            }

          if(P[n].Type == 0)
            {
              P[n]          = P[NumGas - 1];
              P[NumGas - 1] = P[NumPart - 1];
#ifdef BLACKHOLES
              if(P[NumPart-1].Type==5)
                BPP(NumPart-1).PID = NumGas - 1;
#endif
#ifdef STARS
              if(P[NumPart-1].Type==4)
                SPP(NumPart-1).PID = NumGas - 1;
#endif
              Key[n]          = Key[NumGas - 1];
              Key[NumGas - 1] = Key[NumPart - 1];

              SphP[n] = SphP[NumGas - 1];

              NumGas--;
            }
#ifdef BLACKHOLES
          else if(P[n].Type == 5)
            {
              BPP(n) = BhP[NumBhs-1];
              PPB(NumBhs-1).BhID = P[n].BhID; 

              if(n == NumPart-1)
                { 
                  NumBhs--;
                  NumPart--;
                  n--;
                  continue; 
                }  

              P[n] = P[NumPart-1];
              if(P[NumPart-1].Type == 5)
                BPP(NumPart-1).PID = n;

              Key[n] = Key[NumPart - 1];
              
              NumBhs--;
            }
#endif
#ifdef STARS
          else if(P[n].Type == 4)
            {
              SPP(n) = SP[NumStars-1];
              PPS(NumStars-1).SID = P[n].SID; 

              if(n == NumPart-1)
                { 
                  NumStars--;
                  NumPart--;
                  n--;
                  continue; 
                }  

              P[n] = P[NumPart-1];
              if(P[NumPart-1].Type == 4)
                SPP(NumPart-1).PID = n;

              Key[n] = Key[NumPart - 1];
              
              NumStars--;
            }
#endif
          else
            {
              P[n]   = P[NumPart - 1];
#ifdef BLACKHOLES
              if(P[NumPart-1].Type == 5)
                BPP(NumPart-1).PID = n;
#endif
#ifdef STARS
              if(P[NumPart-1].Type == 4)
                SPP(NumPart-1).PID = n;
#endif
              Key[n] = Key[NumPart - 1];
            }

          NumPart--;
          n--;

        } /* target != ThisTask */
    }     /* n < NumPart */

  /**** now resize the storage for the P[] and SphP[] arrays if needed ****/
#ifdef BLACKHOLES
  domain_resize_storage_bhs(count_get_bhs);
#endif
#ifdef STARS
  domain_resize_storage_stars(count_get_stars);
#endif
  domain_resize_storage(count_get, count_get_sph, 1);


  /*****  space has been created, now can do the actual exchange *****/
  int count_totget = count_get_sph;

  if(count_totget)
    {
      memmove(P + NumGas + count_totget, P + NumGas, (NumPart - NumGas) * sizeof(struct particle_data));
      memmove(Key + NumGas + count_totget, Key + NumGas, (NumPart - NumGas) * sizeof(peanokey));

#ifdef BLACKHOLES
      for(i=0; i<NumBhs; i++)
        BhP[i].PID+=count_totget;
#endif  

#ifdef STARS
      for(i=0; i<NumStars; i++)
        SP[i].PID+=count_totget;
#endif  
    }

  for(i = 0; i < NTask; i++)
    {
      count_recv_sph[i] = toGetSph[i];
      count_recv[i]     = toGet[i] - toGetSph[i];
#ifdef BLACKHOLES
      count_recv_bhs[i]  = toGetBhs[i];
#endif
#ifdef STARS
      count_recv_stars[i]  = toGetStars[i];
#endif
    }

  int prec_count;
  for(i = 1, offset_recv_sph[0] = NumGas; i < NTask; i++)
    offset_recv_sph[i] = offset_recv_sph[i - 1] + count_recv_sph[i - 1];
#ifdef BLACKHOLES
  for(i = 1, offset_recv_bhs[0] = NumBhs; i < NTask; i++)
    offset_recv_bhs[i] = offset_recv_bhs[i - 1] + count_recv_bhs[i - 1];
#endif
#ifdef STARS
  for(i = 1, offset_recv_stars[0] = NumStars; i < NTask; i++)
    offset_recv_stars[i] = offset_recv_stars[i - 1] + count_recv_stars[i - 1];
#endif
  prec_count = NumGas + count_get_sph;

  offset_recv[0] = NumPart - NumGas + prec_count;

  for(i = 1; i < NTask; i++)
    offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];

#ifndef USE_MPIALLTOALLV_IN_DOMAINDECOMP

  int ngrp;
#ifdef NO_ISEND_IRECV_IN_DOMAIN /* synchronous communication */
  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_sph[target] > 0 || count_recv_sph[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data), MPI_BYTE, target,
                           TAG_PDATA_SPH, P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data), MPI_BYTE,
                           target, TAG_PDATA_SPH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
                           TAG_SPHDATA, SphP + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct sph_particle_data),
                           MPI_BYTE, target, TAG_SPHDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                           Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
#ifdef BLACKHOLES
if(count_bhs[target] > 0 || count_recv_bhs[target] > 0)
            {
              MPI_Sendrecv(bhBuf + offset_bhs[target], count_bhs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target,
                           TAG_BHDATA, BhP + offset_recv_bhs[target], count_recv_bhs[target] * sizeof(struct bh_particle_data),
                           MPI_BYTE, target, TAG_BHDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
#endif
#ifdef STARS
if(count_stars[target] > 0 || count_recv_stars[target] > 0)
            {
              MPI_Sendrecv(sBuf + offset_stars[target], count_stars[target] * sizeof(struct star_particle_data), MPI_BYTE, target,
                           TAG_BHDATA, SP + offset_recv_stars[target], count_recv_stars[target] * sizeof(struct star_particle_data),
                           MPI_BYTE, target, TAG_STARDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
#endif
          if(count[target] > 0 || count_recv[target] > 0)
            {
              MPI_Sendrecv(partBuf + offset[target], count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                           P + offset_recv[target], count_recv[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                           MPI_COMM_WORLD, MPI_STATUS_IGNORE);

              MPI_Sendrecv(keyBuf + offset[target], count[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY,
                           Key + offset_recv[target], count_recv[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD,
                           MPI_STATUS_IGNORE);
            }
        }
    }

#else  /* #ifdef NO_ISEND_IRECV_IN_DOMAIN */
  /* asynchronous communication */

  MPI_Request *requests = (MPI_Request *)mymalloc_movable(&requests, "requests", 30 * NTask * sizeof(MPI_Request));
  int n_requests        = 0;

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_recv_sph[target] > 0)
            {
              MPI_Irecv(P + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct particle_data), MPI_BYTE, target,
                        TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(SphP + offset_recv_sph[target], count_recv_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
                        TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv_sph[target], count_recv_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                        MPI_COMM_WORLD, &requests[n_requests++]);
            }
#ifdef BLACKHOLES
          if(count_recv_bhs[target] > 0)
            {
              MPI_Irecv(BhP + offset_recv_bhs[target], count_recv_bhs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target,
                        TAG_BHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
#ifdef STARS
          if(count_recv_stars[target] > 0)
            {
              MPI_Irecv(SP + offset_recv_stars[target], count_recv_stars[target] * sizeof(struct star_particle_data), MPI_BYTE, target,
                        TAG_STARDATA, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
          if(count_recv[target] > 0)
            {
              MPI_Irecv(P + offset_recv[target], count_recv[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Irecv(Key + offset_recv[target], count_recv[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD,
                        &requests[n_requests++]);
            }
        }
    }

  MPI_Barrier(MPI_COMM_WORLD); /* not really necessary, but this will guarantee that all receives are
                                  posted before the sends, which helps the stability of MPI on
                                  bluegene, and perhaps some mpich1-clusters */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
        {
          if(count_sph[target] > 0)
            {
              MPI_Isend(partBuf + offset_sph[target], count_sph[target] * sizeof(struct particle_data), MPI_BYTE, target,
                        TAG_PDATA_SPH, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(sphBuf + offset_sph[target], count_sph[target] * sizeof(struct sph_particle_data), MPI_BYTE, target,
                        TAG_SPHDATA, MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(keyBuf + offset_sph[target], count_sph[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY_SPH,
                        MPI_COMM_WORLD, &requests[n_requests++]);
            }
#ifdef BLACKHOLES
          if(count_bhs[target] > 0)
            {
              MPI_Isend(bhBuf + offset_bhs[target], count_bhs[target] * sizeof(struct bh_particle_data), MPI_BYTE, target,
                        TAG_BHDATA, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
#ifdef STARS
          if(count_stars[target] > 0)
            {
              MPI_Isend(sBuf + offset_stars[target], count_stars[target] * sizeof(struct star_particle_data), MPI_BYTE, target,
                        TAG_STARDATA, MPI_COMM_WORLD, &requests[n_requests++]);
            }
#endif
          if(count[target] > 0)
            {
              MPI_Isend(partBuf + offset[target], count[target] * sizeof(struct particle_data), MPI_BYTE, target, TAG_PDATA,
                        MPI_COMM_WORLD, &requests[n_requests++]);

              MPI_Isend(keyBuf + offset[target], count[target] * sizeof(peanokey), MPI_BYTE, target, TAG_KEY, MPI_COMM_WORLD,
                        &requests[n_requests++]);
            }
        }
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#endif /* #ifdef NO_ISEND_IRECV_IN_DOMAIN #else */

#else /* #ifndef USE_MPIALLTOALLV_IN_DOMAINDECOMP */
  /* begins block of myMPI_Alltoallv communications */

  myMPI_Alltoallv(partBuf, count_sph, offset_sph, P, count_recv_sph, offset_recv_sph, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(sphBuf, count_sph, offset_sph, SphP, count_recv_sph, offset_recv_sph, sizeof(struct sph_particle_data), 0,
                  MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count_sph, offset_sph, Key, count_recv_sph, offset_recv_sph, sizeof(peanokey), 0, MPI_COMM_WORLD);

#ifdef BLACKHOLES
 myMPI_Alltoallv(bhBuf, count_bhs, offset_bhs, BhP, count_recv_bhs, offset_recv_bhs, sizeof(struct bh_particle_data), 0,
                  MPI_COMM_WORLD);
#endif

#ifdef STARS
 myMPI_Alltoallv(sBuf, count_stars, offset_stars, SP, count_recv_stars, offset_recv_stars, sizeof(struct star_particle_data), 0,
                  MPI_COMM_WORLD);
#endif

  myMPI_Alltoallv(partBuf, count, offset, P, count_recv, offset_recv, sizeof(struct particle_data), 0, MPI_COMM_WORLD);

  myMPI_Alltoallv(keyBuf, count, offset, Key, count_recv, offset_recv, sizeof(peanokey), 0, MPI_COMM_WORLD);

#endif /* #ifndef USE_MPIALLTOALLV_IN_DOMAINDECOMP #else */
       /* close block of myMPI_Alltoallv communications */

#ifdef BLACKHOLES
  for(int i = NumPart + count_get_sph, j=NumBhs; i < NumPart + count_get; i++)
    { 
      if(P[i].Type == 5)
        {
          P[i].BhID = j;  
          BhP[j].PID = i;
          j++;
        } 
    }
  
  NumBhs += count_get_bhs;
#endif

#ifdef STARS
  for(int i = NumPart + count_get_sph, j=NumStars; i < NumPart + count_get; i++)
    { 
      if(P[i].Type == 4)
        {
          P[i].SID = j;  
          SP[j].PID = i;
          j++;
        } 
    }
  
  NumStars += count_get_stars;
#endif

  NumPart += count_get;
  NumGas += count_get_sph;

  myfree(keyBuf);
#ifdef STARS
  myfree(sBuf);
#endif
#ifdef BLACKHOLES
  myfree(bhBuf);
#endif
  myfree(sphBuf);
  myfree(partBuf);
#ifdef STARS
  myfree(offset_recv_stars);
  myfree(count_recv_stars);
  myfree(offset_stars);
  myfree(count_stars);
#endif
#ifdef BLACKHOLES
  myfree(offset_recv_bhs);
  myfree(count_recv_bhs);
  myfree(offset_bhs);
  myfree(count_bhs);
#endif
  myfree(offset_recv_sph);
  myfree(offset_recv);
  myfree(count_recv_sph);
  myfree(count_recv);
  myfree(offset_sph);
  myfree(offset);
  myfree(count_sph);
  myfree(count);

  double t1 = second();
  mpi_printf("DOMAIN: exchange of %lld particles done. (took %g sec)\n", sumtogo, timediff(t0, t1));
}
