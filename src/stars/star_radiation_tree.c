#include <math.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../stars/star_radiation.h"

static inline int ray_box_intersect(double *ray_pos, double *ray_dir, double *center, double len, double *t_enter, double *t_exit)
{
    double half = 0.5 * len;
    double tmin = -MAX_REAL_NUMBER, tmax = MAX_REAL_NUMBER;

    for(int i = 0; i < 3; i++)
    {
        double lo = center[i] - half;
        double hi = center[i] + half;

        if(fabs(ray_dir[i]) < 1e-12) /* ray parallel to this slab */
        {
            if(ray_pos[i] < lo || ray_pos[i] > hi)
                return 0; /* outside slab, no intersection */
        }
        else
        {
            double inv_dir = 1.0 / ray_dir[i];
            double t1 = (lo - ray_pos[i]) * inv_dir;
            double t2 = (hi - ray_pos[i]) * inv_dir;

            if(t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }

            tmin = (t1 > tmin) ? t1 : tmin;
            tmax = (t2 < tmax) ? t2 : tmax;

            if(tmin > tmax) return 0;
        }
    }

    if(tmax < 0) return 0; /* box is behind the ray */

    *t_enter = tmin;
    *t_exit  = tmax;
    return 1;
}

static inline int ray_sphere_intersect(const double dx, const double dy, const double dz, const double *dir, const double r2, double *t_enter, double *t_exit)
{
  /* Check if ray origin is inside the sphere first.
   * dx/dy/dz = sphere_centre - ray_origin */
  double dist2 = dx*dx + dy*dy + dz*dz;
  if(dist2 < r2)
    {
      /* Origin is inside — find the single forward exit point.
       * Project centre onto ray, then offset by half-chord. */
      double t_closest = dx*dir[0] + dy*dir[1] + dz*dir[2];
      double cx = dx - t_closest*dir[0];
      double cy = dy - t_closest*dir[1];
      double cz = dz - t_closest*dir[2];
      double b2 = cx*cx + cy*cy + cz*cz;
      double dt = sqrt(r2 - b2);
      *t_enter = 0.0;         
      *t_exit  = t_closest + dt;
      return 1;
    }

  /* Sphere centre is outside and ahead */
  double t_closest = dx*dir[0] + dy*dir[1] + dz*dir[2];
  if(t_closest <= 0) return 0;

  double cx = dx - t_closest*dir[0];
  double cy = dy - t_closest*dir[1];
  double cz = dz - t_closest*dir[2];
  double b2 = cx*cx + cy*cy + cz*cz;
  if(b2 >= r2) return 0;

  double dt = sqrt(r2 - b2);
  *t_enter = t_closest - dt;
  *t_exit  = t_closest + dt;
  return 1;
}

static inline int ray_absorb(RayPacket *ray, double chord_length, double density, double kappa[WAVEBANDS], double absorbed_RAD[WAVEBANDS], double *dtau_IR)
{
  for(int w = 0; w < WAVEBANDS; w++)
    {
      absorbed_RAD[w] = 0.0;

      if(!(ray->active_bands & (1u << w)))
        continue;

      double dtau = kappa[w] * density * chord_length;
      double absorbed = ray->RAD[w] * (1.0 - exp(-dtau));

      absorbed_RAD[w] += absorbed;
      ray->RAD[w] -= absorbed;  

      /* deactivate band if it has fallen below the dead-fraction threshold */
      if(ray->RAD[w] < RAD_TRUNC_FRAC * ray->RAD_Initial[w])
        ray->active_bands &= (uint8_t)(~(1u << w));
    }
    
    // IR re-absorption tau
    *dtau_IR = kappa[INFRARED] * density * chord_length;

  return ray->active_bands != 0;
}

/* 
Tree indices are organized as follows:

[0 ... Tree_MaxPart-1] -> real particles

[Tree_MaxPart ... Tree_MaxPart+Tree_MaxNodes-1] -> internal nodes

    └── [Tree_MaxPart ... Tree_FirstNonTopLevelNode-1] -> top-level nodes (replicated everywhere)
 
                       └──  [Tree_FirstNonTopLevelNode ... Tree_MaxPart+Tree_MaxNodes-1] -> local branch nodes

[Tree_MaxPart+Tree_MaxNodes ... Tree_MaxPart+Tree_MaxNodes+NTopleaves-1] -> pseudo-particles

[Tree_MaxPart+Tree_MaxNodes+NTopleaves ... ] -> imported points (Tree_ImportedNodeOffset = Tree_MaxPart + Tree_MaxNodes + NTopleaves)
*/

void raytrace_treewalk(RayPacket *ray, int mode, int target_node, RayExportBuffer *export_buf)
{
  /* local stack for ordering within this domain */
  StackEntry stack[RAY_STACK_SIZE];
  int stack_top = 0;

  /* push entry point */
  if(mode == 0)
    stack[stack_top++] = (StackEntry){0.0, MAX_REAL_NUMBER, Tree_MaxPart}; /* root */
  else
    {
      /* restore whatever was pending from previous rank */
      //if(ray->n_pending >= RAY_STACK_SIZE - 1)
      //  terminate("Pending stack too large to restore!");

      memcpy(stack, ray->pending, ray->n_pending * sizeof(StackEntry));
      stack_top = ray->n_pending;
      ray->n_pending = 0;
      /* push the target node on top - it goes first */
      stack[stack_top++] = (StackEntry){ray->t_enter, ray->t_exit, target_node};
    }
  
  while(stack_top > 0)
    {
      StackEntry cur = stack[--stack_top];
      int no = cur.node;

      /* ---- real particle ---- */
      if(no < Tree_MaxPart)
        { 
          //if(P[no].Type != 0) // this shouldn't really happen in the stack
          //  continue; 
              
          //double px = Tree_Pos_list[3 * no + 0] - ray->pos[0];
          //double py = Tree_Pos_list[3 * no + 1] - ray->pos[1];
          //double pz = Tree_Pos_list[3 * no + 2] - ray->pos[2];
              
          /* assume cell is spherical */
          //double r = cbrt((3.0 * SphP[no].Volume) / (4.0 * M_PI));
          //double r2 = r * r;
              
          /* project onto ray for closest approach */
          //double t_enter, t_exit;
          //if(!ray_sphere_intersect(px, py, pz, ray->dir, r2, &t_enter, &t_exit))
          //  terminate("Should not happen!");
              
          double chord_length = cur.t_exit - cur.t_enter;;
          double density = SphP[no].Density;
              
          double kappa[WAVEBANDS];
          for(int w = 0; w < WAVEBANDS; w++)
            kappa[w] = SphP[no].Kappa[w];
          
          double absorbed[WAVEBANDS];

          double dtau_IR;

          int still_alive = ray_absorb(ray, chord_length, density, kappa, absorbed, &dtau_IR);

          /* deposit absorbed energy into cell, one band at a time */

          // Ionizing Photons
          SphP[no].RAD[IONIZING_H_PHOTONS] += absorbed[IONIZING_H_PHOTONS];
          
          // Ionizing Energy
          double dp = 0.0;
          SphP[no].RAD[IONIZING] += absorbed[IONIZING];
          dp += absorbed[IONIZING] / (CLIGHT / All.UnitVelocity_in_cm_per_s);
              
          double dp_rerad = 0.0;
          for(int w = 2; w < WAVEBANDS; w++)
            {
              SphP[no].RAD[w] += absorbed[w];
              dp_rerad += absorbed[w] * (1 + dtau_IR * ReradiatedFraction[w]) / (CLIGHT / All.UnitVelocity_in_cm_per_s);
            }
         
          SphP[no].StarMomentumFeed[0] += (dp + dp_rerad) * ray->dir[0];
          SphP[no].StarMomentumFeed[1] += (dp + dp_rerad) * ray->dir[1];
          SphP[no].StarMomentumFeed[2] += (dp + dp_rerad) * ray->dir[2];

          if(!still_alive) return; /* all bands are exhausted */     
        }
      /* ---- internal node ---- */
      else if(no < Tree_MaxPart + Tree_MaxNodes)
        {
          struct NODE *nop = &Nodes[no];

          //double t_enter, t_exit;
          //if(!ray_box_intersect(ray->pos, ray->dir, nop->center, nop->len, &t_enter, &t_exit))
          //  terminate("Should not happen!");
          
          // Approx. barnes hut like criterion -> might need for expensive sims 
          //double dtau_node = Kappa * nop->u.d.density * nop->len;
          //if(dtau_node < TAU_OPEN_CRITERION)
          //  {
          //  }

          /* enumerate direct children, sort by t_enter, push */
          StackEntry children[8];
          int nchildren = 0;

          int child = nop->u.d.nextnode;
          while(child != nop->u.d.sibling && child >= 0)
            {
              double ct_enter, ct_exit;
              int hit = 0;

              if(child < Tree_MaxPart) /* cell */
                {
                  if(P[child].Type == 0)
                    {
                      double px = Tree_Pos_list[3 * child + 0] - ray->pos[0];
                      double py = Tree_Pos_list[3 * child + 1] - ray->pos[1];
                      double pz = Tree_Pos_list[3 * child + 2] - ray->pos[2];
                  
                      double r  = cbrt((3.0 * SphP[child].Volume) / (4.0 * M_PI));
                      double r2 = r * r;
                      
                      hit = ray_sphere_intersect(px, py, pz, ray->dir, r2, &ct_enter, &ct_exit);
                    }               
                }
              else if(child < Tree_MaxPart + Tree_MaxNodes) /* internal node */
                {
                  if(Nodes[child].u.d.density > 0)
                    hit = ray_box_intersect(ray->pos, ray->dir, Nodes[child].center, Nodes[child].len, &ct_enter, &ct_exit);
                }
              else if(child >= Tree_ImportedNodeOffset) /* imported point */
                {
                  int n = child - Tree_ImportedNodeOffset;

                  if(Tree_Points[n].Type == 0)
                    {    
                      double px = Tree_Points[n].Pos[0] - ray->pos[0];
                      double py = Tree_Points[n].Pos[1] - ray->pos[1];
                      double pz = Tree_Points[n].Pos[2] - ray->pos[2];

                      double r = cbrt((3.0*Tree_Points[n].Volume)/(4.0*M_PI));
                      double r2 = r * r;

                      hit = ray_sphere_intersect(px, py, pz, ray->dir, r2, &ct_enter, &ct_exit);
                    }
                }
              else /* pseudo-particle - requires export */
                {
                  int pseudo_idx = child - (Tree_MaxPart + Tree_MaxNodes);
                  int top_node = DomainNodeIndex[pseudo_idx];

                  if(Nodes[top_node].u.d.density > 0)
                    hit = ray_box_intersect(ray->pos, ray->dir, Nodes[top_node].center, Nodes[top_node].len, &ct_enter, &ct_exit);
                }

              if(hit)
                {
                  if(nchildren >= 8)
                    terminate("Too many children!");

                  children[nchildren++] = (StackEntry){ct_enter, ct_exit, child};
                }

              /* advance to next direct child via sibling */
              if(child < Tree_MaxPart)
                child = Nextnode[child];
              else if(child < Tree_MaxPart + Tree_MaxNodes)
                child = Nodes[child].u.d.sibling;
              else
                child = Nextnode[child - Tree_MaxNodes];
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
      /* ---- imported particle: remote domain ---- */  
      else if(no >= Tree_ImportedNodeOffset) 
        {
          int n = no - Tree_ImportedNodeOffset;
              
          //double px = Tree_Points[n].Pos[0] - ray->pos[0];
          //double py = Tree_Points[n].Pos[1] - ray->pos[1];
          //double pz = Tree_Points[n].Pos[2] - ray->pos[2];

          //double r  = cbrt((3.0*Tree_Points[n].Volume)/(4.0*M_PI));
          //double r2 = r * r;
          
          //double t_enter, t_exit;
          //if(!ray_sphere_intersect(px, py, pz, ray->dir, r2, &t_enter, &t_exit))
          //  terminate("Should not happen!");
              
          double chord_length = cur.t_exit - cur.t_enter;
          double density = Tree_Points[n].Density;
          
          double kappa[WAVEBANDS];
          for(int w = 0; w < WAVEBANDS; w++)
            kappa[w] = Tree_Points[n].Kappa[w];
          
          double absorbed[WAVEBANDS];

          double dtau_IR;

          int still_alive = ray_absorb(ray, chord_length, density, kappa, absorbed, &dtau_IR);

          /* deposit absorbed energy into cell, one band at a time */

          // Ionizing Photons
          Tree_Points[n].RAD[IONIZING_H_PHOTONS] += absorbed[IONIZING_H_PHOTONS];
          
          // Ionizing Energy
          double dp = 0.0;
          Tree_Points[n].RAD[IONIZING] += absorbed[IONIZING];
          dp += absorbed[IONIZING] / (CLIGHT / All.UnitVelocity_in_cm_per_s);
              
          double dp_rerad = 0.0;
          for(int w = 2; w < WAVEBANDS; w++)
            {
              Tree_Points[n].RAD[w] += absorbed[w];
              dp_rerad += absorbed[w] * (1 + dtau_IR * ReradiatedFraction[w]) / (CLIGHT / All.UnitVelocity_in_cm_per_s);
            }
         
          Tree_Points[n].StarMomentumFeed[0] += (dp + dp_rerad) * ray->dir[0];
          Tree_Points[n].StarMomentumFeed[1] += (dp + dp_rerad) * ray->dir[1];
          Tree_Points[n].StarMomentumFeed[2] += (dp + dp_rerad) * ray->dir[2];

          if(!still_alive) return;  
        }
      /* ---- pseudo-particle: remote domain ---- */  
      else
        {
          int task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)];
          int remote_node = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];

          /* pack pending */
          if(stack_top >= RAY_STACK_SIZE - 1) 
            terminate("Too many pending entries to export!");

          ray->n_pending = stack_top;
          memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

          ray->t_enter = cur.t_enter;
          ray->t_exit = cur.t_exit;
          ray->target_node = remote_node;  /* store for mode==1 */

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
    }
}