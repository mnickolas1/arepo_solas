#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"


/*
 * Exact ray transport on the Voronoi mesh via face connectivity.
 *
 * For a ray x(t) = x0 + t n inside cell i, cell membership is the
 * intersection of the bisector half spaces against every face-defining
 * neighbour j in DC(i):
 *
 *     (x - m_ij) . d_ij <= 0 ,  d_ij = s_j - s_i ,  m_ij = (s_i + s_j)/2
 *
 * Substituting x(t) and solving for equality,
 *
 *     t_j = [ (m_ij - x0) . d_ij ] / ( n . d_ij ) ,     n . d_ij > 0
 *
 * The exit is t_exit = min_j t_j over forward-facing bisectors only, and the
 * argmin IS the next cell. Backward-facing bisectors have n . d_ij < 0 and
 * drop out automatically - in particular the face we just entered through,
 * whose sign flips exactly (not approximately) on the crossing.
 *
 * All quantities are held relative to the current generator.
 */

/* Sphere-equivalent cross section of a cell, pi r^2, used for the ray
   splitting criterion. With this convention All.RaySplitFactor reads as the
   minimum number of rays required to sample each cell */
#define RAY_CELL_CROSS_SECTION(r) (M_PI * (r) * (r))

static inline double cell_dtau(int i, double length, double N_H2_ray,
                               ChannelsDtau dtau[WAVEBANDS])
{
  const double dust_length = SphP[i].OpacityScaling[CH_DUST] * length;
  
  const double dN_H2 = SphP[i].OpacityScaling[CH_H2] * length;  
  const double dtau_line = h2shield_dtau(N_H2_ray, dN_H2);

  double ionizing_length[3];
  for(int s = 0; s < 3; s++)
    ionizing_length[s] = SphP[i].OpacityScaling[CH_HI + s] * length;

  for(int w = 0; w < WAVEBANDS; w++)
    {
      for(int c = 0; c < CHANNELS; c++)
        dtau[w].E[c] = dtau[w].N[c] = 0.0;

      const uint8_t ch = BandChannels[w];

      if(ch & (1u << CH_DUST))
        {
          dtau[w].E[CH_DUST] = Kappa_E[w] * dust_length;
          dtau[w].N[CH_DUST] = Kappa_N[w] * dust_length;
        }

      if(ch & (1u << CH_H2))
        dtau[w].E[CH_H2] = dtau[w].N[CH_H2] = dtau_line;

      for(int s = 0; s < 3; s++)
        if(ch & (1u << (CH_HI + s)))
          {
            dtau[w].E[CH_HI + s] = Sigma_E[w - IONIZING_HI][s] * ionizing_length[s];
            dtau[w].N[CH_HI + s] = Sigma_N[w - IONIZING_HI][s] * ionizing_length[s];
          }
    }

  return dN_H2;
}

static inline int ray_absorb(RayPacket *ray, const ChannelsDtau dtau[WAVEBANDS],
                             Absorption *a)
{
  memset(a, 0, sizeof(*a));

  for(int w = 0; w < WAVEBANDS; w++)
    {
      if(ray->Radiated[w].Energy < RAD_TRUNC_FRAC * ray->Radiated_Init[w].Energy &&
         ray->Radiated[w].Photons < RAD_TRUNC_FRAC * ray->Radiated_Init[w].Photons)
        ray->active_bands &= (uint8_t)(~(1u << w));

      if(!(ray->active_bands & (1u << w)))
        continue;

      double tot_E = 0.0, tot_N = 0.0;
      for(int c = 0; c < CHANNELS; c++)
        {
          tot_E += dtau[w].E[c];
          tot_N += dtau[w].N[c];
        }

      if(tot_E <= 0.0 && tot_N <= 0.0)
        continue;

      /* expm1, not 1-exp: more accurate for optically thin cells */
      double dE = ray->Radiated[w].Energy * -expm1(-tot_E);
      double dN = ray->Radiated[w].Photons * -expm1(-tot_N);

      a->Band[w].Energy = dE;
      a->Band[w].Photons = dN;

      ray->Radiated[w].Energy -= dE;
      ray->Radiated[w].Photons -= dN;

      for(int c = 0; c < CHANNELS; c++)
        {
          if(tot_E > 0.0) a->Ch[w][c].Energy = dE * dtau[w].E[c] / tot_E;
          if(tot_N > 0.0) a->Ch[w][c].Photons = dN * dtau[w].N[c] / tot_N;
        }
    }

  return ray->active_bands != 0;
}

/* Deposition into one cell over an exact path length */
/* Returns 0 once every band of this ray is exhausted */
static inline int ray_deposit(RayPacket *ray, int i, double length)
{
  if(length <= 0.0)
    return ray->active_bands != 0;

  ChannelsDtau dtau[WAVEBANDS];
  double dN_H2 = cell_dtau(i, length, ray->N_H2, dtau);

  /* Process ray */
  Absorption a;
  int still_alive = ray_absorb(ray, dtau, &a);

  /* Accumulate H2 column */
  ray->N_H2 += dN_H2;

  /* Reradiation in the IR (boosts momentum) */
  double Dtau_IR = dtau_IR(i, length);

  const double c_code = CLIGHT / All.cf_UnitVelocity_in_cm_per_s;

  const double mj = P[i].Mass + SphP[i].StarMassFeed;

  double mom[3], smf[3];
  for(int k = 0; k < 3; k++)
    {
      mom[k] = SphP[i].Momentum[k];
      smf[k] = SphP[i].StarMomentumFeed[k];
    }

  /* Deposit absorbed energy into the cell, one band at a time */
  double dK_total = 0.0;

  for(int w = 0; w < WAVEBANDS; w++)
    {
      const double E_w = a.Band[w].Energy;
      if(E_w <= 0.0)
        continue;

      /* Only the dust channel reradiates in the IR */
      const double f_dust = a.Ch[w][CH_DUST].Energy / E_w;

      double dp = E_w * (1.0 + f_dust * Dtau_IR * ReradiatedFraction[w]) / c_code / All.cf_atime;

      double dp_vec[3] = {dp * ray->dir[0], dp * ray->dir[1], dp * ray->dir[2]};

      /* Partially updated state */
      double pj[3];
      for(int k = 0; k < 3; k++)
        pj[k] = mom[k] + smf[k];

      double sq_momentum = dp_vec[0] * dp_vec[0] + dp_vec[1] * dp_vec[1] + dp_vec[2] * dp_vec[2];

      double cross = 2 * (pj[0] * dp_vec[0] + pj[1] * dp_vec[1] + pj[2] * dp_vec[2]);

      double dK = (sq_momentum + cross) / (2 * mj);

      smf[0] += dp_vec[0];
      smf[1] += dp_vec[1];
      smf[2] += dp_vec[2];

      const double f_K = dK / E_w;
      for(int c = 0; c < CHANNELS; c++)
        a.Ch[w][c].Energy *= (1.0 - f_K);

      dK_total += dK;
    }

#ifdef RADIATION_PRESSURE
  for(int k = 0; k < 3; k++)
    SphP[i].StarMomentumFeed[k] = smf[k];

  SphP[i].StarEnergyFeed += dK_total;
  All.StarFeedbackLocal[2] += dK_total;
#endif

#ifdef PHOTOELECTRIC_HEATING
  SphP[i].AbsorbedPE += a.Ch[ULTRAVIOLET][CH_DUST].Energy * TrueAbsorbedFraction[ULTRAVIOLET]
                      + a.Ch[LYMAN_WERNER][CH_DUST].Energy * TrueAbsorbedFraction[LYMAN_WERNER];
#endif

#ifdef DISSOCIATION
  SphP[i].AbsorbedH2Line += a.Ch[LYMAN_WERNER][CH_H2].Photons;
#endif

#ifdef PHOTOIONIZATION
  for(int s = 0; s < 3; s++)
    {
      for(int w = IONIZING_HI; w <= IONIZING_HeII; w++)
        {
          SphP[i].AbsorbedIonizing[s].Energy += a.Ch[w][CH_HI + s].Energy;
          SphP[i].AbsorbedIonizing[s].Photons += a.Ch[w][CH_HI + s].Photons;
        }
    }
#endif

  return still_alive;
}

/* Non-periodic wall: the face-defining "neighbour" is the 
   mirror image of this cell across a box face 
   Both REFLECTIVE_* = 1 and = 2 build the same mirror ghost; 
   for now we treat either as outflow and drop the ray */
static inline int dc_is_boundary(int q)
{
#if defined(REFLECTIVE_X) || defined(REFLECTIVE_Y) || defined(REFLECTIVE_Z)
  const int flags = DC[q].image_flags & MASK;

#ifdef REFLECTIVE_X
  if(flags & (MASK_X_SHIFT_RIGHT | MASK_X_SHIFT_LEFT))
    return 1;
#endif

#ifdef REFLECTIVE_Y
  if(flags & (MASK_Y_SHIFT_RIGHT | MASK_Y_SHIFT_LEFT))
    return 1;
#endif

#ifdef REFLECTIVE_Z
  if(flags & (MASK_Z_SHIFT_RIGHT | MASK_Z_SHIFT_LEFT))
    return 1;
#endif
#endif

  return 0;
}

#define RAY_MAX_LOCATE 4096

/* Returns 0 if the head now lies inside ray->cell (possibly on its boundary),
   1 if the ray was handed to another rank and the caller must return */
static int voronoi_relocate(RayPacket *ray, RayComms *comm)
{
  for(int it = 0; it < RAY_MAX_LOCATE; it++)
    {
      int q_best = -1;
      double v_best = 0.0;
      double d_best[3] = {0.0, 0.0, 0.0};

      const int i = ray->cell;

      const double r_cell = get_cell_radius(i);
      const double eps = RAY_TOL * r_cell;

      const double sx = P[i].Pos[0], sy = P[i].Pos[1], sz = P[i].Pos[2];
      
      const double px = ray->pos[0] + ray->t * ray->dir[0];
      const double py = ray->pos[1] + ray->t * ray->dir[1];
      const double pz = ray->pos[2] + ray->t * ray->dir[2];

      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          int dp = DC[q].dp_index;
          int particle = Mesh.DP[dp].index;
          
          if(particle >= 0)
            {
              double dx = Mesh.DP[dp].x - sx;
              double dy = Mesh.DP[dp].y - sy;
              double dz = Mesh.DP[dp].z - sz;

              double d2 = dx*dx + dy*dy + dz*dz;

              /* v = 1/2 (|x-s_i|^2 - |x-s_j|^2); 
                 v/|d| is the signed distance to the bisector plane 
                 so that v > 0 means s_j is the closer generator */
              double v = (px*dx + py*dy + pz*dz) - 0.5 * d2;

              if(v > v_best && v > eps * sqrt(d2))
                {
                  v_best = v;
                  q_best = q;
                  
                  d_best[0] = dx; 
                  d_best[1] = dy; 
                  d_best[2] = dz;
                }
            }

          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
        }

      if(q_best < 0)
        {
          ray->locate_head = 0;
          return 0;                  
        }

      /* Head is outside the box: the child has already escaped */
      if(dc_is_boundary(q_best))
        return 1;

      ray->pos[0] -= d_best[0];
      ray->pos[1] -= d_best[1];
      ray->pos[2] -= d_best[2];

      ray->cell = DC[q_best].index;

      if(DC[q_best].task != ThisTask)
        {
          append_export(comm, ray, DC[q_best].task);
          return 1;                   
        }
    }

  terminate("voronoi_relocate(): locate walk did not converge for ray!\n");
  return 0;
}

/* Exit face search */
/* P. Camps (2013)
 * Min-reduction of t_j over the forward-facing bisectors of cell i
 * Returns the DC index of the exit connection (or -1 if no forward-facing bisector exists, i.e. the cell is unbounded along n)
 */
static inline int voronoi_exit_face(const RayPacket *ray, int i, double eps, double *t_out, double d_out[3])
{
  int q_best = -1;
  double t_best = MAX_REAL_NUMBER;
  double area_best = -1.0;
  int area_best_valid = 0;

  const double nx = ray->dir[0], ny = ray->dir[1], nz = ray->dir[2];
  const double sx = P[i].Pos[0], sy = P[i].Pos[1], sz = P[i].Pos[2];

  const double px = ray->pos[0] + ray->t * nx;
  const double py = ray->pos[1] + ray->t * ny;
  const double pz = ray->pos[2] + ray->t * nz;

  int q = SphP[i].first_connection;

  while(q >= 0)
    {
      if(q >= MaxNvc)
        terminate("voronoi_exit_face(): strange connectivity q=%d MaxNvc=%d cell=%d!\n", q, MaxNvc, i);

      int dp = DC[q].dp_index;
      int particle = Mesh.DP[dp].index;
          
      /* Cell has been removed */
      if(particle < 0)
        {
          if(q == SphP[i].last_connection)
            break;

          q = DC[q].next;
          
          continue;
        }

      /* Generator offset; carries the periodic/reflective image shift */
      double dx = Mesh.DP[dp].x - sx;
      double dy = Mesh.DP[dp].y - sy;
      double dz = Mesh.DP[dp].z - sz;

      double ndotd = nx * dx + ny * dy + nz * dz;

      /* Backward-facing bisector: non-binding, includes the entry face */
      if(ndotd > 0.0)
        {
          /* (m - x0) . d with m = d/2 in this frame */
          double num = 0.5 * (dx*dx + dy*dy + dz*dz) - (px * dx + py * dy + pz * dz);

          double t = num / ndotd;

          /* Genuinely closer */
          if(t < t_best - eps)
            {
              t_best = t;
              q_best = q;

              d_out[0] = dx;
              d_out[1] = dy;
              d_out[2] = dz;

              area_best_valid = 0;
            }
          /* Within error */
          else if(t < t_best + eps && q_best >= 0)
            {
              if(!area_best_valid)
                {
                  area_best = Mesh.VF[DC[q_best].vf_index].area;
                  area_best_valid = 1;
                }

              double area = Mesh.VF[DC[q].vf_index].area;

              /* Tie break */
              if(area > area_best)
                {
                  if(t < t_best)
                    t_best = t;

                  q_best = q;
                  area_best = area;

                  d_out[0] = dx;
                  d_out[1] = dy;
                  d_out[2] = dz;
                }
            }
        }

      if(q == SphP[i].last_connection)
        break;

      q = DC[q].next;
    }

  *t_out = t_best;

  return q_best;
}

/* Main traversal */
void raytrace_voronoi(RayPacket *ray, RayWorkStack *work, RayComms *comm)
{
  long long steps = 0;

  while(1)
    {
      if(ray->locate_head)
        {
          if(voronoi_relocate(ray, comm))
            return;
        }

      int i = ray->cell;

      const double r_cell = get_cell_radius(i);
      const double eps = RAY_TOL * r_cell;

      /* Adaptive splitting, evaluated on cell entry */
      if(ray->nside < NSIDE_MAX && ray->t > 0.0)
        {
          double A_cell = RAY_CELL_CROSS_SECTION(r_cell);
          double Omega_ray = 4.0 * M_PI / (12.0 * (double)ray->nside * (double)ray->nside);

          if(A_cell / (ray->t * ray->t) < All.RaySplitFactor * Omega_ray)
            {
              RayPacket children[4];
              split_ray(ray, children);
              
              for(int k = 0; k < 4; k++)
                {
                  if(voronoi_relocate(&children[k], comm))
                    continue;           
    
                  append_ray(work, &children[k]);
                }
              
              return;
            }
        }

      /* Exact path length through this cell */
      double t_step, d[3];
      int q = voronoi_exit_face(ray, i, eps, &t_step, d);

      if(q < 0)
        {
          terminate("raytrace_voronoi(): cell %d unbounded along ray - stale DC list!\n", i);
        }

      if(t_step < 0.0)
        {
          if(t_step < -eps)
            warn("raytrace_voronoi(): cell %d gives negative t %g along ray - stale DC list?\n", i, t_step);
          
          t_step = 0.0;
        }

      int truncated = 0;
      if(ray->t + t_step >= ray->t_maximum)
        {
          t_step = ray->t_maximum - ray->t;
          truncated = 1;
        }

      /* Absorption, heating, radiation pressure */
      int still_alive = ray_deposit(ray, i, t_step);

      SphP[i].RTCost += 1.0f;

      ray->t += t_step;

      if(!still_alive || truncated)
        return;

      /* Outflow boundary */
      if(dc_is_boundary(q))
        return;

      /* Re-anchor on the neighbour's generator */
      ray->pos[0] -= d[0];
      ray->pos[1] -= d[1];
      ray->pos[2] -= d[2];

      int task = DC[q].task;

      ray->cell = DC[q].index;

      /* Hand off to the owning rank */
      if(task != ThisTask)
        {
          append_export(comm, ray, task);
          return;
        }

      /* A self-connection means the exit face is a reflective mirror image of itself (not treated here) */
      if(ray->cell == i)
        terminate("raytrace_voronoi(): self-connection at cell %d (mirror boundary)!\n", i);

      if(++steps > RAY_MAX_CELL_STEPS)
        {       
          warn("raytrace_voronoi(): ray exceeded %d cell steps on task %d?\n", RAY_MAX_CELL_STEPS, ThisTask);          
          return;
        }

      if((steps & RAY_STEPS_PROGRESS_MASK) == 0)
        ray_comms_progress(comm);
    }
}