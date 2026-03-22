#ifndef NGBTREE_H
#define NGBTREE_H

/*! Variables for neighbor tree */
extern int Ngb_MaxPart;
extern int Ngb_NumNodes;
extern int Ngb_MaxNodes;
extern int Ngb_FirstNonTopLevelNode;
extern int Ngb_NextFreeNode;
extern int *Ngb_Father;
extern int *Ngb_Marker;
extern int Ngb_MarkerValue;

extern int *Ngb_DomainNodeIndex;
extern int *DomainListOfLocalTopleaves;
extern int *DomainNLocalTopleave;
extern int *DomainFirstLocTopleave;
extern int *Ngb_Nextnode;

/*! The ngb-tree data structure
 */
extern struct NgbNODE
{
  union
  {
    int suns[8]; /*!< temporary pointers to daughter nodes */
    struct
    {
      int sibling;
      int nextnode;
      MyNgbTreeFloat range_min[3];
      MyNgbTreeFloat range_max[3];
    } d;
  } u;

  MyNgbTreeFloat vertex_vmin[3];
  MyNgbTreeFloat vertex_vmax[3];

  int father;

  integertime Ti_Current;
} * Ngb_Nodes;

extern struct ExtNgbNODE
{
  float vmin[3];
  float vmax[3];
  float MaxCsnd;
} * ExtNgb_Nodes;

#ifdef STAR_RADIATION_ACTIVE
extern struct RtNgbNODE
{
    /* geometric cell bounds — opening angle + intersection */
    float center[3];
    float len;

    /* RT quantities — volume-weighted, accumulated bottom-up */
    int nchildren;

    float volume;
    float density_kappa[WAVEBANDS];
} * RtNgb_Nodes;
#endif

#endif