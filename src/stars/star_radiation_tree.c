#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"

static inline int ray_box_intersect(double *ray_pos, double *ray_dir, MyNgbTreeFloat *rmin, MyNgbTreeFloat *rmax, double *t_enter, double *t_exit)
{
  double tmin = -MAX_REAL_NUMBER, tmax = MAX_REAL_NUMBER;

  for(int i = 0; i < 3; i++)
    {
      if(fabs(ray_dir[i]) < 1e-12) /* ray parallel to this slab */
        {
          if(ray_pos[i] < rmin[i] || ray_pos[i] > rmax[i])
            return 0;
        }
      else
        {
          double inv_dir = 1.0 / ray_dir[i];
          double t1 = (rmin[i] - ray_pos[i]) * inv_dir;
          double t2 = (rmax[i] - ray_pos[i]) * inv_dir;

          if(t1 > t2) 
            { 
              double tmp = t1; 
              t1 = t2; 
              t2 = tmp; 
            }

          tmin = t1 > tmin ? t1 : tmin;
          tmax = t2 < tmax ? t2 : tmax;

          if(tmin > tmax) 
            return 0;
        }
    }

    if(tmax < 0) 
      return 0;

    *t_enter = fmax(tmin, 0.0);
    *t_exit  = tmax;

    return 1;
}

static inline int ray_sphere_intersect(const double dx, const double dy, const double dz, const double *dir, const double r2, double *t_enter, double *t_exit)
{
  /* Check if ray origin is inside the sphere first.
   * dx/dy/dz = sphere_centre - ray_origin */
  double dist2 = dx * dx + dy * dy + dz * dz;
  
  if(dist2 < r2)
    {
      /* Origin is inside — find the single forward exit point.
       * Project centre onto ray, then offset by half-chord. */
      double t_closest = dx * dir[0] + dy * dir[1] + dz * dir[2];
      double cx = dx - t_closest * dir[0];
      double cy = dy - t_closest * dir[1];
      double cz = dz - t_closest * dir[2];
      double b2 = cx * cx + cy * cy + cz * cz;
      double dt = sqrt(r2 - b2);
      *t_enter = 0.0;         
      *t_exit  = t_closest + dt;
      return 1;
    }

  /* Sphere centre is outside and ahead */
  double t_closest = dx * dir[0] + dy * dir[1] + dz * dir[2];
  
  if(t_closest <= 0) 
    return 0;

  double cx = dx - t_closest * dir[0];
  double cy = dy - t_closest * dir[1];
  double cz = dz - t_closest * dir[2];
  double b2 = cx * cx + cy * cy + cz * cz;
  
  if(b2 >= r2) 
    return 0;

  double dt = sqrt(r2 - b2);
  *t_enter = t_closest - dt;
  *t_exit  = t_closest + dt;
  
  return 1;
}

static inline int ray_absorb(RayPacket *ray, double chord_length, double density_kappa[WAVEBANDS], double absorbed_RAD[WAVEBANDS], double *dtau_IR)
{
  for(int w = 0; w < WAVEBANDS; w++)
    {
      absorbed_RAD[w] = 0.0;

      if(!(ray->active_bands & (1u << w)))
        continue;

      double dtau = density_kappa[w] * chord_length;
      double absorbed = ray->RAD[w] * (1.0 - exp(-dtau));

      absorbed_RAD[w] += absorbed;
      ray->RAD[w] -= absorbed;  

      /* deactivate band if it has fallen below the dead-fraction threshold */
      if(ray->RAD[w] < RAD_TRUNC_FRAC * ray->RAD_Initial[w])
        ray->active_bands &= (uint8_t)(~(1u << w));
    }
    
    // IR re-absorption tau
    *dtau_IR = density_kappa[INFRARED] * chord_length;

  return ray->active_bands != 0;
}

/* 
Tree indices are organized as follows:

[0 ... Ngb_MaxPart-1] -> real particles

[Ngb_MaxPart ... Ngb_MaxPart+Ngb_MaxNodes-1] -> internal nodes

    └── [Ngb_MaxPart ... Tree_FirstNonTopLevelNode-1] -> top-level nodes (replicated everywhere)
 
                       └──  [Tree_FirstNonTopLevelNode ... Ngb_MaxPart+Ngb_MaxNodes-1] -> local branch nodes

[Ngb_MaxPart+Ngb_MaxNodes ... Ngb_MaxPart+Ngb_MaxNodes+NTopleaves-1] -> pseudo-particles
*/

void raytrace_treewalk(RayPacket *ray, RayWorkStack *work, RayExportBuffer *export_buf)
{
  double xtmp, ytmp, ztmp;
  
  /* local stack for ordering within this domain */
  StackEntry stack[RAY_STACK_SIZE];
  int stack_top = 0;

  /* push entry point */
  if(ray->target_node < 0 )
    stack[stack_top++] = (StackEntry){0.0, MAX_REAL_NUMBER, Ngb_MaxPart}; /* root */
  else
    {
      memcpy(stack, ray->pending, ray->n_pending * sizeof(StackEntry));
      stack_top = ray->n_pending;
      ray->n_pending = 0;
      /* push the target node on top - it goes first */
      stack[stack_top++] = (StackEntry){ray->t, ray->t_exit, ray->target_node};
    }
  
  while(stack_top > 0)
    {
      StackEntry cur = stack[--stack_top];
      int no = cur.node;

      /* ---- real particle ---- */
      if(no < Ngb_MaxPart)
        {     
          double chord_length = cur.t_exit - cur.t_enter;
              
          double density_kappa[WAVEBANDS];
          for(int w = 0; w < WAVEBANDS; w++)
            density_kappa[w] = SphP[no].Density * SphP[no].Kappa[w];
          
          double absorbed[WAVEBANDS];
          double dtau_IR;

          int still_alive = ray_absorb(ray, chord_length, density_kappa, absorbed, &dtau_IR);

          /* deposit absorbed energy into cell, one band at a time */

          // Ionizing Photons
          SphP[no].RAD[IONIZING_H_PHOTONS] += absorbed[IONIZING_H_PHOTONS];
          
          // Ionizing Energy
          double dp = 0.0;
          SphP[no].RAD[IONIZING] += absorbed[IONIZING];
          dp += absorbed[IONIZING] / (CLIGHT / All.cf_UnitVelocity_in_cm_per_s) / All.cf_atime;
              
          double dp_rerad = 0.0;
          for(int w = 2; w < WAVEBANDS; w++)
            {
              SphP[no].RAD[w] += absorbed[w];
              dp_rerad += absorbed[w] * (1 + dtau_IR * ReradiatedFraction[w]) / (CLIGHT / All.cf_UnitVelocity_in_cm_per_s) / All.cf_atime;
            }
         
          SphP[no].StarMomentumFeed[0] += (dp + dp_rerad) * ray->dir[0];
          SphP[no].StarMomentumFeed[1] += (dp + dp_rerad) * ray->dir[1];
          SphP[no].StarMomentumFeed[2] += (dp + dp_rerad) * ray->dir[2];

          ray->t = cur.t_exit;

          if(ray->t == ray->t_maximum) 
            {
              ray->is_paused = 1; 
              return;
            }
          
          /* all bands are exhausted */
          if(!still_alive) 
            return;      
        }
      /* ---- internal node ---- */
      else if(no < Ngb_MaxPart + Ngb_MaxNodes)
        {
          struct NgbNODE *nop = &Ngb_Nodes[no];

#ifdef RAD_OPENING_ANGLE
          /* This should only trigger for non-top level nodes */ 
          if (no >= Ngb_FirstNonTopLevelNode)
            {
              /* --- Barnes-Hut opening criterion --- */
              double cx = 0.5 * (nop->u.d.range_max[0] + nop->u.d.range_min[0]);
              double cy = 0.5 * (nop->u.d.range_max[1] + nop->u.d.range_min[1]);
              double cz = 0.5 * (nop->u.d.range_max[2] + nop->u.d.range_min[2]);

              double ddx = cx - ray->pos[0];
              double ddy = cy - ray->pos[1];
              double ddz = cz - ray->pos[2];
              double dist2 = ddx*ddx + ddy*ddy + ddz*ddz;

              double dx = nop->u.d.range_max[0] - nop->u.d.range_min[0];
              double dy = nop->u.d.range_max[1] - nop->u.d.range_min[1];
              double dz = nop->u.d.range_max[2] - nop->u.d.range_min[2];
              double len2 = dx*dx + dy*dy + dz*dz;
          
              /* node is far enough — treat as single slab */
              if(dist2 > 0 && len2 / dist2 < All.RadOpeningAngle * All.RadOpeningAngle)
                {
                  double chord_length = cur.t_exit - cur.t_enter;
              
                  double density_kappa[WAVEBANDS];
                  for(int w = 0; w < WAVEBANDS; w++)
                    density_kappa[w] = RtNgb_Nodes[no].density_kappa[w];  /* volume-weighted mean kappa * density */

                  double absorbed[WAVEBANDS];
                  double dtau_IR;

                  int still_alive = ray_absorb(ray, chord_length, density_kappa, absorbed, &dtau_IR);

                  /* accumulate for later distribution to children */
                  for(int w = 0; w < WAVEBANDS; w++)
                    RtNgb_Nodes[no].RAD[w] += absorbed[w];

                  ray->t = cur.t_exit;

                  if(ray->t == ray->t_maximum) 
                    {
                      ray->is_paused = 1; 
                      return;
                    }

                  if(!still_alive) 
                    return;
                
                  continue;  /* don't open this node */
                }
            }
#endif      
          /* Adaptive splitting criterion */
          if(ray->nside < NSIDE_MAX)
            {
              double cx = 0.5 * (nop->u.d.range_max[0] + nop->u.d.range_min[0]);
              double cy = 0.5 * (nop->u.d.range_max[1] + nop->u.d.range_min[1]);
              double cz = 0.5 * (nop->u.d.range_max[2] + nop->u.d.range_min[2]);

              double ddx = cx - ray->pos[0];
              double ddy = cy - ray->pos[1];
              double ddz = cz - ray->pos[2];
              double dist2 = ddx*ddx + ddy*ddy + ddz*ddz;

              double dx = nop->u.d.range_max[0] - nop->u.d.range_min[0];
              double dy = nop->u.d.range_max[1] - nop->u.d.range_min[1];
              double dz = nop->u.d.range_max[2] - nop->u.d.range_min[2];
              double len2 = dx*dx + dy*dy + dz*dz;

              /* use number of actual children for adaptive f */
             
              double f_eff = fmax(1.0, (double)RtNgb_Nodes[no].nchildren); /* at least 1 ray per child */

              double coeff = 3.0 * f_eff / M_PI; /* from Ω_ray = 4π/(12·nside²) */

              if(dist2 > 0.0 && len2 * coeff * (double)(ray->nside * ray->nside) > dist2)
                {
                  /* ray is too coarse — push split children to split_buf, consume parent */
                  RayPacket children[4];
                  if(split_ray(ray, children))
                    {
                      /* pack pending */
                      if(stack_top >= RAY_STACK_SIZE - 1) 
                        terminate("Too many pending entries to split!");

                      ray->n_pending = stack_top;
                      memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

                      ray->t = cur.t_enter;
                      ray->t_exit = cur.t_exit;
                      ray->target_node = no;  /* store for re-entry */

                      for(int k = 0; k < 4; k++)
                        {
                          if(work->n >= work->capacity)
                            terminate("Work stack overflow!");
                          
                          append_ray(work, &children[k]);
                        }
                      /* parent ray is consumed */
                      return;   
                    }
                  /* else: at NSIDE_MAX, fall through and open the node normally */
                }
            }
          /* open node and enumerate children -> sort by t_enter, push */
          StackEntry children[8];
          int nchildren = 0;

          int child = nop->u.d.nextnode;
          while(child != nop->u.d.sibling && child >= 0)
            {
              double t_enter, t_exit;
              int hit = 0;

              if(child < Ngb_MaxPart) /* cell */
                {
                  double px = NEAREST_X(P[child].Pos[0] - ray->pos[0]);
                  double py = NEAREST_Y(P[child].Pos[1] - ray->pos[1]);
                  double pz = NEAREST_Z(P[child].Pos[2] - ray->pos[2]);
                  
                  double r  = cbrt((3.0 * SphP[child].Volume) / (4.0 * M_PI));
                  double r2 = r * r;
                      
                  hit = ray_sphere_intersect(px, py, pz, ray->dir, r2, &t_enter, &t_exit);              
                }
              else if(child < Ngb_MaxPart + Ngb_MaxNodes) /* internal node */
                hit = ray_box_intersect(ray->pos, ray->dir, Ngb_Nodes[child].u.d.range_min, Ngb_Nodes[child].u.d.range_max, &t_enter, &t_exit);
            
              else /* pseudo-particle - requires export */
                {
                  int pseudo_idx = child - (Ngb_MaxPart + Ngb_MaxNodes);
                  int top_node = Ngb_DomainNodeIndex[pseudo_idx];

                  hit = ray_box_intersect(ray->pos, ray->dir, Ngb_Nodes[top_node].u.d.range_min, Ngb_Nodes[top_node].u.d.range_max, &t_enter, &t_exit);
                }

              if(hit)
                {
                  if(t_enter > ray->t_maximum)
                  /* child is beyond max travel distance - skip entirely */
                    continue;
                  else
                    {
                      t_exit = fmin(t_exit, ray->t_maximum);  /* limit traversal distance */
                  
                      if(nchildren >= 8)
                        terminate("Too many children!");

                      children[nchildren++] = (StackEntry){t_enter, t_exit, child};
                    }
                }

              /* advance to next direct child via sibling */
              if(child < Ngb_MaxPart)
                child = Ngb_Nextnode[child];
              else if(child < Ngb_MaxPart + Ngb_MaxNodes)
                child = Ngb_Nodes[child].u.d.sibling;
              else
                child = Ngb_Nextnode[child - Ngb_MaxNodes];
            }

          /* sort ascending by t_enter */
          for(int i = 1; i < nchildren; i++)
            {
              StackEntry key = children[i];
              int j = i - 1;
              while(j >= 0 && children[j].t_enter > key.t_enter)
                {
                  children[j+1] = children[j];
                  j--;
                }
              children[j+1] = key;
            }

          /* push in reverse so smallest t_enter is popped first */
          for(int i = nchildren - 1; i >= 0; i--)
            {
              if(stack_top >= RAY_STACK_SIZE)
              terminate("Ray stack overflow!");

              stack[stack_top++] = children[i];
            }
        }
    
      /* ---- pseudo-particle: remote domain ---- */  
      else
        {
          int task = DomainTask[no - (Ngb_MaxPart + Ngb_MaxNodes)];
          int remote_node = Ngb_DomainNodeIndex[no - (Ngb_MaxPart + Ngb_MaxNodes)];

          /* pack pending */
          if(stack_top >= RAY_STACK_SIZE - 1) 
            terminate("Too many pending entries to export!");

          ray->n_pending = stack_top;
          memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

          ray->t = cur.t_enter;
          ray->t_exit = cur.t_exit;
          ray->target_node = remote_node;  /* store for re-entry */

          /* add to export buffer */
          if(export_buf->n < export_buf->capacity)
            {
              export_buf->rays[export_buf->n] = *ray;
              export_buf->task[export_buf->n] = task;
              export_buf->n++;
            }
          else
            terminate("Export buffer full!");

        /* this rank is done with this ray */
        return;
        }
       
      if(stack_top == 0)
        if(ray->t < ray->t_maximum) 
          {
            ray->t = ray->t_maximum;
            ray->is_paused = 1; 
            return;
          }
    }
}