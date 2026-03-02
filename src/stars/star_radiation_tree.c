
/* Wavebands */
typedef enum
{
  LYMAN_WERNER = 0,
  ULTRAVIOLET,
  OPTICAL,
  INFRARED,
  WB_COUNT
} Waveband;


/* Opacity coefficients */
static double Kappa[WB_COUNT] = {
  1.0,  // LW
  1.0,  // UV
  1.0,  // OP
  1.0   // IR
};

#define RAY_STACK_SIZE 64

typedef struct 
{
  double t_enter;
  int node;
} StackEntry;

typedef struct 
{
  double pos[3];
  double dir[3];
  double RAD_Ionizing;
  double t;
  int ray_id;
  int home_task;
  /* pending top-level nodes still to traverse after current domain */
  StackEntry pending[RAY_STACK_SIZE];
  int  n_pending;
} RayData;

static int ray_box_intersect(double *ray_pos, double *ray_dir, double *center, double len, double *t_enter, double *t_exit)
{
    double half = 0.5 * len;
    double tmin = -1e30, tmax = 1e30;

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

/* 
Tree indices are saved as follows:

[0 ... Tree_MaxPart-1] -> real particles

[Tree_MaxPart ... Tree_MaxPart+Tree_MaxNodes-1] -> internal nodes

    └── [Tree_MaxPart ... Tree_FirstNonTopLevelNode-1] -> top-level nodes (replicated everywhere)
 
                       └──  [Tree_FirstNonTopLevelNode ... Tree_MaxPart+Tree_MaxNodes-1] -> local branch nodes

[Tree_MaxPart+Tree_MaxNodes ... Tree_MaxPart+Tree_MaxNodes+NTopleaves-1] -> pseudo-particles

[Tree_MaxPart+Tree_MaxNodes+NTopleaves ... ] -> imported points (Tree_ImportedNodeOffset = Tree_MaxPart + Tree_MaxNodes + NTopleaves)
*/

void raytrace_treewalk(RayData *ray, int mode, int target_node)
{
  /* local stack for ordering within this domain */
  StackEntry stack[RAY_STACK_SIZE];
  int stack_top = 0;

  /* push entry point */
  if(mode == 0)
    stack[stack_top++] = (StackEntry){0.0, Tree_MaxPart}; /* root */
  else
    {
      /* restore whatever was pending from previous rank */
      memcpy(stack, ray->pending, ray->n_pending * sizeof(StackEntry));
      stack_top = ray->n_pending;
      ray->n_pending = 0;
      /* push the target node on top - it goes first */
      stack[stack_top++] = (StackEntry){ray->t, target_node};
    }
  
  while(stack_top > 0)
    {
      StackEntry cur = stack[--stack_top];
      int no = cur.node;

      /* ---- real particle ---- */
      if(no < Tree_MaxPart)
        {
          double px = Tree_Pos_list[3*no+0] - ray->pos[0];
          double py = Tree_Pos_list[3*no+1] - ray->pos[1];
          double pz = Tree_Pos_list[3*no+2] - ray->pos[2];
                  
          /* assume cell is a sphere */
          /* project onto ray for closest approach */
          double t_closest = px*ray->dir[0] + py*ray->dir[1] + pz*ray->dir[2];

          /* perpendicular distance squared from ray to sphere center */
          double cx = px - t_closest*ray->dir[0];
          double cy = py - t_closest*ray->dir[1];
          double cz = pz - t_closest*ray->dir[2];
          double b2 = cx*cx + cy*cy + cz*cz;

          double r  = cbrt((3.0*SphP[no].Volume)/(4.0*M_PI)); 
          double r2 = r * r;

          if(b2 < r2 && t_closest > 0)
            {
              /* ray intersects sphere */
              double chord_length = 2.0 * sqrt(r2 - b2); /* path through sphere */

              double density = SphP[no].Density;
              double metallicity = SphP[no].Metallicity;
              double kappa = Kappa[LYMAN_WERNER]; 

              /* optical depth */
              double dtau = kappa * density * metallicity * chord_length;

              /* absorption */
              SphP[no].RAD_Ionizing += ray->RAD_Ionizing * (1 - exp(-dtau));

              /* attenuation */
              ray->RAD_Ionizing *= exp(-dtau);

              /* early termination */
              if(ray->RAD_Ionizing < RAD_BACKGROUND)
                return;
            }
          else
            terminate("Should not happen!");
        }
      /* ---- internal node ---- */
      else if(no < Tree_MaxPart + Tree_MaxNodes)
        {
          /* mode==1 guard: don't escape back into top-level tree */
          if(mode == 1 && no < Tree_FirstNonTopLevelNode)
            terminate("Should not happen!");

          struct NODE *nop = &Nodes[no];

          double t_enter, t_exit;
          if(!ray_box_intersect(ray->pos, ray->dir, nop->center, nop->len, &t_enter, &t_exit))
            terminate("Should not happen!");

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
                  double px = Tree_Pos_list[3*child+0] - ray->pos[0];
                  double py = Tree_Pos_list[3*child+1] - ray->pos[1];
                  double pz = Tree_Pos_list[3*child+2] - ray->pos[2];
                  
                  /* assume cell is a sphere */
                  /* project onto ray for closest approach */
                  double t_closest = px*ray->dir[0] + py*ray->dir[1] + pz*ray->dir[2];

                  /* perpendicular distance squared from ray to sphere center */
                  double cx = px - t_closest*ray->dir[0];
                  double cy = py - t_closest*ray->dir[1];
                  double cz = pz - t_closest*ray->dir[2];
                  double b2 = cx*cx + cy*cy + cz*cz;

                  double r  = cbrt((3.0*SphP[child].Volume)/(4.0*M_PI)); 
                  double r2 = r * r;

                  if(b2 < r2 && t_closest > 0)
                    {
                      /* ray intersects sphere */
                      double dt = sqrt(r2 - b2); /* half chord length */
                      ct_enter = t_closest - dt;
                      ct_exit = t_closest + dt;
                      hit = 1;
                    }
                  else
                    hit = 0;
                }
              else if(child < Tree_MaxPart + Tree_MaxNodes) /* internal node */
                {
                  hit = ray_box_intersect(ray->pos, ray->dir, Nodes[child].center, Nodes[child].len, &ct_enter, &ct_exit);
                }
              else if(child >= Tree_ImportedNodeOffset) /* imported point */
                {
                  int n = child - Tree_ImportedNodeOffset;
    
                  double px = Tree_Points[n].Pos[0] - ray->pos[0];
                  double py = Tree_Points[n].Pos[1] - ray->pos[1];
                  double pz = Tree_Points[n].Pos[2] - ray->pos[2];

                  double t_closest = px*ray->dir[0] + py*ray->dir[1] + pz*ray->dir[2];

                  double cx = px - t_closest*ray->dir[0];
                  double cy = py - t_closest*ray->dir[1];
                  double cz = pz - t_closest*ray->dir[2];
                  double b2 = cx*cx + cy*cy + cz*cz;

                  double r = cbrt((3.0*Tree_Points[n].Volume)/(4.0*M_PI));
                  double r2 = r * r;

                  if(b2 < r2 && t_closest > 0)
                    {
                      /* ray intersects sphere */
                      double dt = sqrt(r2 - b2); /* half chord length */
                      ct_enter = t_closest - dt;
                      ct_exit = t_closest + dt;
                      hit = 1;
                    }
                }
              else /* pseudo-particle - requires export */
                {
                  int pseudo_idx = child - (Tree_MaxPart + Tree_MaxNodes);
                  int top_node = DomainNodeIndex[pseudo_idx];

                  hit = ray_box_intersect(ray->pos, ray->dir, Nodes[top_node].center, Nodes[top_node].len, &ct_enter, &ct_exit);
                }

              if(hit)
                children[nchildren++] = (StackEntry){ct_enter, child};

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
            stack[stack_top++] = children[i];
        }

      /* ---- imported particle: remote domain ---- */  
      else if(no >= Tree_ImportedNodeOffset) 
        {
          int n = no - Tree_ImportedNodeOffset;
    
          double px = Tree_Points[n].Pos[0] - ray->pos[0];
          double py = Tree_Points[n].Pos[1] - ray->pos[1];
          double pz = Tree_Points[n].Pos[2] - ray->pos[2];
        
          /* assume cell is a sphere */
          /* project onto ray for closest approach */
          double t_closest = px*ray->dir[0] + py*ray->dir[1] + pz*ray->dir[2];

          /* perpendicular distance squared from ray to sphere center */
          double cx = px - t_closest*ray->dir[0];
          double cy = py - t_closest*ray->dir[1];
          double cz = pz - t_closest*ray->dir[2];
          double b2 = cx*cx + cy*cy + cz*cz;

          double r  = cbrt((3.0*Tree_Points[n].Volume)/(4.0*M_PI));
          double r2 = r * r;

          if(b2 < r2 && t_closest > 0)
            {
              /* ray intersects sphere */
              double chord_length = 2.0 * sqrt(r2 - b2); /* path through sphere */

              double density = Tree_Points[n].Density;
              double metallicity = Tree_Points[n].Metallicity;
              double kappa = Kappa[LYMAN_WERNER]; 

              /* optical depth */
              double dtau = kappa * density * metallicity * chord_length;

              /* absorption */
              Tree_Points[n].RAD_Ionizing += ray->RAD_Ionizing * (1 - exp(-dtau));

              /* attenuation */
              ray->RAD_Ionizing *= exp(-dtau);

              /* early termination */
              if(ray->RAD_Ionizing < RAD_BACKGROUND)
                return;
            }
          else
            terminate("Should not happen!");
        }
      /* ---- pseudo-particle: remote domain ---- */  
      else
        {
          int task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)];
          int remote_node = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];

          /* everything left on local stack is downstream of this domain,
             but we need to check: are any remaining stack entries 
             actually local? No - because we sorted by t_enter and 
             this pseudo-particle came before them, so they are all 
             further along the ray = further domains */

          /* pack pending: current local stack + ray's existing pending */
          ray->n_pending = stack_top;
          memcpy(ray->pending, stack, stack_top * sizeof(StackEntry));

          /* forward ray to remote rank */
          ray->t = cur.t_enter;
          export_ray(ray, task, remote_node);

          /* this rank is done with this ray */
          return;
        }
    }

  /* ray truly finished - send result home */
  send_result_home(ray);
}

/*! \brief Prepares node to be exported.
 *
 *  \param[in] no Index of node.
 *  \param[in] i Index of particle.
 *  \param[in] thread_id ID of thread.
 *
 *  \return 0
 */
int tree_treefind_export_node_threads(int no, int i, int thread_id)
{
  /* The task indicated by the pseudoparticle node */
  int task = DomainNewTask[no - (Tree_MaxPart + Tree_MaxNodes)];

  if(Thread[thread_id].Exportflag[task] != i)
    {
      Thread[thread_id].Exportflag[task]     = i;
      int nexp                               = Thread[thread_id].Nexport++;
      Thread[thread_id].PartList[nexp].Task  = task;
      Thread[thread_id].PartList[nexp].Index = i;
      Thread[thread_id].ExportSpace -= Thread[thread_id].ItemSize;
    }

  int nexp                      = Thread[thread_id].NexportNodes++;
  nexp                          = -1 - nexp;
  struct datanodelist *nodelist = (struct datanodelist *)(((char *)Thread[thread_id].PartList) + Thread[thread_id].InitialSpace);
  nodelist[nexp].Task           = task;
  nodelist[nexp].Index          = i;
  nodelist[nexp].Node           = DomainNodeIndex[no - (Tree_MaxPart + Tree_MaxNodes)];
  Thread[thread_id].ExportSpace -= sizeof(struct datanodelist) + sizeof(int);
  return 0;
}