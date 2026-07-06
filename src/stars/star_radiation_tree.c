#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"


static inline int ray_box_intersect(double *ray_pos, double *ray_dir, MyNgbTreeFloat *rmin, MyNgbTreeFloat *rmax, double *t_enter, double *t_exit)
{
  double tmin = -MAX_REAL_NUMBER, tmax = MAX_REAL_NUMBER;

  for(int i = 0; i < 3; i++)
    {
      /* Ray parallel to this slab */
      if(fabs(ray_dir[i]) < 1e-12) 
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
  /* Check if ray origin is inside the sphere first */
  double dist2 = dx * dx + dy * dy + dz * dz;
  
  if(dist2 < r2)
    {
      /* Origin is inside - find the single forward exit point */
      /* Project centre onto ray, then offset by half-chord */
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

static inline int ray_absorb(RayPacket *ray, double chord_length, double density_kappa[WAVEBANDS], WavebandData absorbed[WAVEBANDS], double *dtau_IR)
{
  for(int w = 0; w < WAVEBANDS; w++)
    {
      absorbed[w].Energy = absorbed[w].Photons = 0.0;

      if(!(ray->active_bands & (1u << w)))
        continue;

      double dtau = density_kappa[w] * chord_length;
      
      double absorbed_energy = ray->Radiated[w].Energy * (1.0 - exp(-dtau));
      double absorbed_photons = ray->Radiated[w].Photons * (1.0 - exp(-dtau));

      absorbed[w].Energy += absorbed_energy;
      ray->Radiated[w].Energy -= absorbed_energy;
      
      absorbed[w].Photons += absorbed_photons;
      ray->Radiated[w].Photons -= absorbed_photons; 

      /* Deactivate band if it has fallen below the dead-fraction threshold */
      if(ray->Radiated[w].Energy < RAD_TRUNC_FRAC * ray->Radiated_Init[w].Energy && 
      ray->Radiated[w].Photons < RAD_TRUNC_FRAC * ray->Radiated_Init[w].Photons)
        ray->active_bands &= (uint8_t)(~(1u << w));
    }
    
    /* IR re-absorption tau */
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
  
  /* Local stack for ordering within this domain */
  StackEntry stack[RAY_STACK_SIZE];
  int stack_top = 0;

  /* Push entry point */
  if(ray->target_node < 0 )
    /* Root */
    stack[stack_top++] = (StackEntry){0.0, MAX_REAL_NUMBER, Ngb_MaxPart}; 
  else
    {
      memcpy(stack, ray->pending, ray->n_pending * sizeof(StackEntry));
      stack_top = ray->n_pending;
      ray->n_pending = 0;
      /* Push the target node on top - it goes first */
      stack[stack_top++] = (StackEntry){ray->t, ray->t_exit, ray->target_node};
    }
  
  while(stack_top > 0)
    {
      StackEntry cur = stack[--stack_top];
      int no = cur.node;

      /* ---- Cell ---- */
      if(no < Ngb_MaxPart)
        {     
          if(P[no].Type != 0 || P[no].Mass == 0 || P[no].ID == 0)
            continue;

          if(P[no].Ti_Current != All.Ti_Current)
            drift_particle(no, All.Ti_Current);

          double chord_length = cur.t_exit - cur.t_enter;
              
          double density_kappa[WAVEBANDS];
          for(int w = 0; w < WAVEBANDS; w++)
            density_kappa[w] = SphP[no].Density * SphP[no].Kappa[w];
          
          WavebandData absorbed[WAVEBANDS];
          double dtau_IR;

          int still_alive = ray_absorb(ray, chord_length, density_kappa, absorbed, &dtau_IR);
          
          /* Deposit absorbed energy into cells, one band at a time */

          double dK_total = 0.0;
          for(int w = 0; w < WAVEBANDS; w++)
            {
              double dp;

              if(w == LYMAN_WERNER || w == IONIZING_HI || w == IONIZING_HeI || w == IONIZING_HeII)
                dp = absorbed[w].Energy / (CLIGHT / All.cf_UnitVelocity_in_cm_per_s) / All.cf_atime;

              else
                dp = absorbed[w].Energy * (1.0 + dtau_IR * ReradiatedFraction[w]) / (CLIGHT / All.cf_UnitVelocity_in_cm_per_s) / All.cf_atime;        
            
              double dp_vec[3] = {dp * ray->dir[0], dp * ray->dir[1], dp * ray->dir[2]};

              /* Partially updated state */
              double mj, pj[3];

              mj = P[no].Mass + SphP[no].StarMassFeed;
              for(int k = 0; k < 3; k++)
                pj[k] = SphP[no].Momentum[k] + SphP[no].StarMomentumFeed[k];

              double sq_momentum = dp_vec[0]*dp_vec[0] + dp_vec[1]*dp_vec[1] + dp_vec[2]*dp_vec[2];

              double cross = 2 * (pj[0] * dp_vec[0] + pj[1] * dp_vec[1] + pj[2] * dp_vec[2]);

              double dK = (sq_momentum + cross) / (2 * mj);

              SphP[no].StarMomentumFeed[0] += dp_vec[0];
              SphP[no].StarMomentumFeed[1] += dp_vec[1];
              SphP[no].StarMomentumFeed[2] += dp_vec[2];

              SphP[no].Absorbed[w].Energy += absorbed[w].Energy - dK;

              dK_total += dK;
            }

          SphP[no].StarEnergyFeed += dK_total;
          All.StarFeedbackLocal[2] += dK_total;
          
          /* Deposit absorbed photons into cells, one band at a time */
          
          /* Dissociating Photons */
          SphP[no].Absorbed[LYMAN_WERNER].Photons += absorbed[LYMAN_WERNER].Photons;

          /* Ionizing Photons */
          SphP[no].Absorbed[IONIZING_HI].Photons += absorbed[IONIZING_HI].Photons;
          SphP[no].Absorbed[IONIZING_HeI].Photons += absorbed[IONIZING_HeI].Photons;
          SphP[no].Absorbed[IONIZING_HeII].Photons += absorbed[IONIZING_HeII].Photons;
         
          ray->t = cur.t_exit;

          if(ray->t == ray->t_maximum) 
            {
              ray->is_paused = 1; 
              return;
            }
          
          /* All bands are exhausted */
          if(!still_alive) 
            return;      
        }
      /* ---- Internal node ---- */
      else if(no < Ngb_MaxPart + Ngb_MaxNodes)
        {
          struct NgbNODE *nop = &Ngb_Nodes[no];
     
          if(nop->Ti_Current != All.Ti_Current)
            drift_node(nop, All.Ti_Current);

#ifdef RAD_OPENING_ANGLE
          /* This should only trigger for non-top level nodes */ 
          if (no >= Ngb_FirstNonTopLevelNode)
            {
              /* -- Barnes-Hut opening criterion -- */
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
          
              /* Node is far enough - treat as single slab */
              if(dist2 > 0 && len2 / dist2 < All.RadOpeningAngle * All.RadOpeningAngle)
                {
                  double chord_length = cur.t_exit - cur.t_enter;
              
                  double density_kappa[WAVEBANDS];
                  for(int w = 0; w < WAVEBANDS; w++)
                    /* Volume-weighted mean kappa * density */
                    density_kappa[w] = RtNgb_Nodes[no].density_kappa[w];  

                  WavebandData absorbed[WAVEBANDS];
                  double dtau_IR;

                  int still_alive = ray_absorb(ray, chord_length, density_kappa, absorbed, &dtau_IR);

                  /* Accumulate for later distribution to children */
                  for(int w = 0; w < WAVEBANDS; w++)
                    {
                      RtNgb_Nodes[no].Absorbed[w].Energy += absorbed[w].Energy;
                      RtNgb_Nodes[no].Absorbed[w].Photons += absorbed[w].Photons;
                    }

                  ray->t = cur.t_exit;

                  if(ray->t == ray->t_maximum) 
                    {
                      ray->is_paused = 1; 
                      return;
                    }

                  if(!still_alive) 
                    return;
                  
                  /* Don't open this node */
                  continue;  
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

              /* Use number of actual children for adaptive f */
             
              /* At least 1 ray per child */
              double f_eff = fmax(1.0, (double)RtNgb_Nodes[no].nchildren); 
              
              /* Omega_ray = 4pi/(12*nside^2) */
              double coeff = 3.0 * f_eff / M_PI; 

              if(dist2 > 0.0 && len2 * coeff * (double)(ray->nside * ray->nside) > dist2)
                {
                  /* Ray is too coarse - push split children to split_buf, consume parent */
                  RayPacket children[4];
                  if(split_ray(ray, children))
                    {
                      /* Pack pending */
                      if(stack_top >= RAY_STACK_SIZE - 1) 
                        terminate("Too many pending entries to split!");

                      ray->n_pending = stack_top;
                      memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

                      ray->t = cur.t_enter;
                      ray->t_exit = cur.t_exit;
                      /* Store for re-entry */
                      ray->target_node = no;  

                      for(int k = 0; k < 4; k++)
                        {
                          if(work->n >= work->capacity)
                            terminate("Work stack overflow!");
                          
                          append_ray(work, &children[k]);
                        }
                      /* Parent ray is consumed */
                      return;   
                    }
                  /* Else: at NSIDE_MAX, fall through and open the node normally */
                }
            }
          /* Open node and enumerate children -> sort by t_enter, push */
          StackEntry children[8];
          int nchildren = 0;

          int child = nop->u.d.nextnode;
          while(child != nop->u.d.sibling && child >= 0)
            {
              double t_enter, t_exit;
              int hit = 0;

              /* Cell */
              if(child < Ngb_MaxPart) 
                {
                  if(P[child].Ti_Current != All.Ti_Current)
                    drift_particle(child, All.Ti_Current);

                  double px = NEAREST_X(P[child].Pos[0] - ray->pos[0]);
                  double py = NEAREST_Y(P[child].Pos[1] - ray->pos[1]);
                  double pz = NEAREST_Z(P[child].Pos[2] - ray->pos[2]);
                  
                  double r = get_cell_radius(child);
                  double r2 = r * r;
                      
                  hit = ray_sphere_intersect(px, py, pz, ray->dir, r2, &t_enter, &t_exit);              
                }
              /* Internal node */  
              else if(child < Ngb_MaxPart + Ngb_MaxNodes) 
                {   
                  struct NgbNODE *child_nop = &Ngb_Nodes[child];

                  if(child_nop->Ti_Current != All.Ti_Current)
                    drift_node(child_nop, All.Ti_Current);

                  hit = ray_box_intersect(ray->pos, ray->dir, child_nop->u.d.range_min, child_nop->u.d.range_max, &t_enter, &t_exit);
                }
              /* Pseudo-particle: remote domain */
              else 
                {
                  int pseudo_idx = child - (Ngb_MaxPart + Ngb_MaxNodes);
                  int top_node = Ngb_DomainNodeIndex[pseudo_idx];

                  struct NgbNODE *pseudo_nop = &Ngb_Nodes[top_node];

                  if(pseudo_nop->Ti_Current != All.Ti_Current)
                    drift_node(pseudo_nop, All.Ti_Current);

                  hit = ray_box_intersect(ray->pos, ray->dir, pseudo_nop->u.d.range_min, pseudo_nop->u.d.range_max, &t_enter, &t_exit);
                }

              if(hit)
                {
                  if(t_enter > ray->t_maximum)
                  /* Child is beyond max travel distance - skip entirely */
                    continue;
                  else
                    {
                      /* Limit traversal distance */
                      t_exit = fmin(t_exit, ray->t_maximum);  
                  
                      if(nchildren >= 8)
                        terminate("Too many children!");

                      children[nchildren++] = (StackEntry){t_enter, t_exit, child};
                    }
                }

              /* Advance to next direct child via sibling */
              if(child < Ngb_MaxPart)
                child = Ngb_Nextnode[child];
              else if(child < Ngb_MaxPart + Ngb_MaxNodes)
                child = Ngb_Nodes[child].u.d.sibling;
              else
                child = Ngb_Nextnode[child - Ngb_MaxNodes];
            }

          /* Sort ascending by t_enter */
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

          /* Push in reverse so smallest t_enter is popped first */
          for(int i = nchildren - 1; i >= 0; i--)
            {
              if(stack_top >= RAY_STACK_SIZE)
              terminate("Ray stack overflow!");

              stack[stack_top++] = children[i];
            }
        }
    
      /* ---- Pseudo-particle: remote domain ---- */  
      else
        {
          int task = DomainTask[no - (Ngb_MaxPart + Ngb_MaxNodes)];
          int remote_node = Ngb_DomainNodeIndex[no - (Ngb_MaxPart + Ngb_MaxNodes)];

          /* Pack pending */
          if(stack_top >= RAY_STACK_SIZE - 1) 
            terminate("Too many pending entries to export!");

          ray->n_pending = stack_top;
          memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

          ray->t = cur.t_enter;
          ray->t_exit = cur.t_exit;
          /* Store for re-entry */
          ray->target_node = remote_node;  

          /* Add to export buffer */
          if(export_buf->n < export_buf->capacity)
            {
              export_buf->rays[export_buf->n] = *ray;
              export_buf->task[export_buf->n] = task;
              export_buf->n++;
            }
          else
            terminate("Export buffer full!");

        /* This rank is done with this ray */
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