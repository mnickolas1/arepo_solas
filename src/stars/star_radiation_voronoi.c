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
 * whose sign flips exactly (not approximately) on the crossing. There is
 * therefore no need to nudge the ray off the face at re-entry.
 *
 * All quantities are held relative to the current generator s_i, so d_ij is
 * literally Mesh.DP[dp] - P[i].Pos and inherits AREPO's image shifts for
 * free. No NEAREST_* wrapping appears below.
 */

/* Sphere-equivalent cross section of a cell, pi r^2, used for the ray
   splitting criterion. With this convention All.RaySplitFactor reads as the
   minimum number of rays required to sample each cell */
#define RAY_CELL_CROSS_SECTION(r) (M_PI * (r) * (r))

static inline int ray_absorb(RayPacket *ray, double Dtau_E[WAVEBANDS], double Dtau_N[WAVEBANDS], WavebandData absorbed[WAVEBANDS],
                             double dN_H2, double *lw_line)
{
  absorbed[IONIZING_HI].Energy = absorbed[IONIZING_HI].Photons = 0.0;

  /* Deactivate band if it has fallen below the dead-fraction threshold */
  if(ray->Radiated[IONIZING_HI].Energy < RAD_TRUNC_FRAC * ray->Radiated_Init[IONIZING_HI].Energy && 
     ray->Radiated[IONIZING_HI].Photons < RAD_TRUNC_FRAC * ray->Radiated_Init[IONIZING_HI].Photons)
     ray->active_bands &= (uint8_t)(~(1u << IONIZING_HI));

    if(!(ray->active_bands & (1u << IONIZING_HI)))
      return ray->active_bands != 0;

    double absorbed_energy = ray->Radiated[IONIZING_HI].Energy * (1.0 - exp(-Dtau_E[IONIZING_HI]));
    double absorbed_photons = ray->Radiated[IONIZING_HI].Photons * (1.0 - exp(-Dtau_N[IONIZING_HI]));

    absorbed[IONIZING_HI].Energy += absorbed_energy;
    ray->Radiated[IONIZING_HI].Energy -= absorbed_energy;
      
    absorbed[IONIZING_HI].Photons += absorbed_photons;
    ray->Radiated[IONIZING_HI].Photons -= absorbed_photons; 

    return ray->active_bands != 0;
}

/* ------------------------------------------------------------------ */
/* Deposition into one cell over an exact path length                  */
/* ------------------------------------------------------------------ */

/* Returns 0 once every band of this ray is exhausted */
static inline int ray_deposit(RayPacket *ray, int i, double length)
{
  if(length <= 0.0)
    return ray->active_bands != 0;

  double Dtau_E[WAVEBANDS], Dtau_N[WAVEBANDS];

  for(int w = 0; w < WAVEBANDS; w++)
    {
      Dtau_E[w] = SphP[i].DtauOverLength_E[w] * length;
      Dtau_N[w] = SphP[i].DtauOverLength_N[w] * length;
    }

  WavebandData absorbed[WAVEBANDS];

  /* Line dissociation */
  double dN_H2 = 0.0;
  /* Fraction of LW absorption that goes into H2 line dissociation */
  double lw_line[2] = {0.0, 0.0};

  /* Process ray */
  int still_alive = ray_absorb(ray, Dtau_E, Dtau_N, absorbed, dN_H2, lw_line);

  SphP[i].Absorbed[IONIZING_HI].Energy += absorbed[IONIZING_HI].Energy;
  SphP[i].Absorbed[IONIZING_HI].Photons += absorbed[IONIZING_HI].Photons;

  return still_alive;
}

#define RAY_MAX_LOCATE 4096

/* Returns 0 if the head now lies inside ray->cell (possibly on its boundary),
   1 if the ray was handed to another rank and the caller must return. */
static int voronoi_relocate(RayPacket *ray, RayExportBuffer *export_buf)
{
  for(int it = 0; it < RAY_MAX_LOCATE; it++)
    {
      int q_best = -1;
      double v_best = 0.0;
      double d_best[3] = {0.0, 0.0, 0.0};

      const int i = ray->cell;

      const double eps = RAY_TOL * get_cell_radius(i);

      const double sx = P[i].Pos[0], sy = P[i].Pos[1], sz = P[i].Pos[2];
      
      const double px = ray->pos[0] + ray->t * ray->dir[0];
      const double py = ray->pos[1] + ray->t * ray->dir[1];
      const double pz = ray->pos[2] + ray->t * ray->dir[2];

      int q = SphP[i].first_connection;
      while(q >= 0)
        {
          int dp = DC[q].dp_index;

          if(Mesh.DP[dp].index >= 0)
            {
              double dx = Mesh.DP[dp].x - sx;
              double dy = Mesh.DP[dp].y - sy;
              double dz = Mesh.DP[dp].z - sz;

              double d2 = dx*dx + dy*dy + dz*dz;

              /* v = 1/2 (|x-s_i|^2 - |x-s_j|^2); v/|d| is the signed distance
                 to the bisector plane. v > 0 means s_j is the closer generator. */
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

      ray->pos[0] -= d_best[0];
      ray->pos[1] -= d_best[1];
      ray->pos[2] -= d_best[2];

      ray->cell = DC[q_best].index;

      if(DC[q_best].task != ThisTask)
        {
          append_export(export_buf, ray, DC[q_best].task);
          return 1;                   
        }
    }

  terminate("RAYTRACE_VORONOI: locate walk did not converge for ray %d (star %llu)\n",
            ray->ray_id, (unsigned long long)ray->star_id);
  return 0;
}

/* ------------------------------------------------------------------ */
/* Exit face search                                                    */
/* ------------------------------------------------------------------ */
/*
 * Min-reduction of t_j over the forward-facing bisectors of cell i.
 * Returns the DC index of the exit connection (or -1 if no forward-facing
 * bisector exists, i.e. the cell is unbounded along n), writes the path
 * length to the exit into *t_out, and writes the winning generator offset
 * d = s_j - s_i into d_out so the caller does not have to re-fetch DP.
 */
static inline int voronoi_exit_face(const RayPacket *ray, int i, double r_cell, double *t_out, double d_out[3])
{
  int q_best = -1;
  double t_best = MAX_REAL_NUMBER;
  double area_best = -1.0;
  int area_best_valid = 0;

  const double eps = RAY_TOL * r_cell;

  const double nx = ray->dir[0], ny = ray->dir[1], nz = ray->dir[2];
  const double sx = P[i].Pos[0], sy = P[i].Pos[1], sz = P[i].Pos[2];

  const double px = ray->pos[0] + ray->t * nx;
  const double py = ray->pos[1] + ray->t * ny;
  const double pz = ray->pos[2] + ray->t * nz;

  int q = SphP[i].first_connection;

  while(q >= 0)
    {
      if(q >= MaxNvc)
        terminate("RAYTRACE_VORONOI: strange connectivity q=%d MaxNvc=%d cell=%d\n", q, MaxNvc, i);

      int dp = DC[q].dp_index;

      /* Cell has been removed */
      if(Mesh.DP[dp].index < 0)
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

          if(t < t_best - eps)
            {
              t_best = t;
              q_best = q;

              d_out[0] = dx;
              d_out[1] = dy;
              d_out[2] = dz;

              /* Defer the VF load until a tie actually needs it */
              area_best_valid = 0;
            }
          else if(t < t_best + eps && q_best >= 0)
            {
              if(!area_best_valid)
                {
                  area_best = Mesh.VF[DC[q_best].vf_index].area;
                  area_best_valid = 1;
                }

              double area = Mesh.VF[DC[q].vf_index].area;

              if(area > area_best)
                {
                  /* Keep the tighter bound, step through the real face */
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

/* ------------------------------------------------------------------ */
/* Main traversal                                                      */
/* ------------------------------------------------------------------ */

void raytrace_voronoi(RayPacket *ray, RayWorkStack *work, RayExportBuffer *export_buf)
{
  long long steps = 0;

  while(1)
    {
      if(ray->locate_head)
        {
          if(voronoi_relocate(ray, export_buf))
            return;
        }

      int i = ray->cell;

      const double r_cell = get_cell_radius(i);

      /* ---- Adaptive splitting, evaluated on cell entry ----
       * Omega_cell / Omega_ray is the number of rays crossing this cell;
       * split while that would fall below All.RaySplitFactor. */
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
                  if(voronoi_relocate(&children[k], export_buf))
                    continue;           
    
                  append_ray(work, &children[k]);
                }
              
              return;
            }
        }

      /* ---- Exact path length through this cell ---- */
      double t_step, d[3];
      int q = voronoi_exit_face(ray, i, r_cell, &t_step, d);

      if(q < 0)
        {
          terminate("RAYTRACE_VORONOI: cell %d unbounded along ray %d (star %llu) - stale DC list?\n", i, ray->ray_id,
                   (unsigned long long)ray->star_id);
        }

      if(t_step < 0.0)
        {
          terminate("RAYTRACE_VORONOI: cell %d gives negative t along ray %d (star %llu) - stale DC list?\n", i, ray->ray_id,
                   (unsigned long long)ray->star_id);
        }

      int truncated = 0;
      if(ray->t + t_step >= ray->t_maximum)
        {
          t_step = ray->t_maximum - ray->t;
          truncated = 1;
        }

      /* ---- Absorption, heating, dissociation, radiation pressure ---- */
      int still_alive = ray_deposit(ray, i, t_step);

      ray->t += t_step;

      if(!still_alive || truncated)
        return;

      /* ---- Re-anchor on the neighbour's generator ----
       * x_new - s_j = (x_old - s_i) + t n - (s_j - s_i)
       * Exact, image-shift aware, and leaves the ray sitting precisely on
       * the shared bisector so the return face is excluded by sign alone. */
      ray->pos[0] -= d[0];
      ray->pos[1] -= d[1];
      ray->pos[2] -= d[2];

      int task = DC[q].task;

      ray->cell = DC[q].index;

      /* Hand off to the owning rank */
      if(task != ThisTask)
        {
          append_export(export_buf, ray, task);
          return;
        }

      /* A self-connection means the exit face is a reflective mirror image of
         this cell, which this scheme does not transport across */
      if(ray->cell == i)
        terminate("RAYTRACE_VORONOI: self-connection at cell %d (mirror boundary?)\n", i);

      if(++steps > RAY_MAX_CELL_STEPS)
        {
          warn("RAYTRACE_VORONOI: ray %d from star %llu exceeded %d cell steps on task %d, dropping\n", ray->ray_id,
               (unsigned long long)ray->star_id, RAY_MAX_CELL_STEPS, ThisTask);
          return;
        }
    }
}