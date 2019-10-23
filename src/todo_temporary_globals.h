#ifndef TODO_TEMPORARY_GLOBALS_H_
#define TODO_TEMPORARY_GLOBALS_H_

#define MLADENASSN 200 /*assume max number of neighbours*/



#include "hydro.h"

/* ================================ */
struct gizmo_debug_dump {
/* ================================ */



  float pos[3];                     /* particle positions */
  long long id;                     /* particle ID */
  float h;                          /* smoothing length */

  float wgrads_store[3];            /* store sum of individual cartesian gradients contributions */
  float volume_store;               /* particle volume */
  float omega;                      /* normalization for psi */

  int nneigh;                       /* number of neighbours this particle interacts with */
  int neighbour_ids[MLADENASSN];    /* IDs of each neighbour this particle interacts with */
  float Aij[2*MLADENASSN];          /* effective surface towards each neighbour */

  int nneigh_grads;
  int neighbour_ids_grad[MLADENASSN];      /* IDs of neighbour for individual gradient contributions */
  float grads_sum_contrib[2*MLADENASSN];   /* contributions to the gradient sum from each neighbour */
  float dwdr[MLADENASSN];                  /* radial derivative of the kernel */

  float grads_sum_dx[2*MLADENASSN];        /* pi.x - pj.x*/
  float r[MLADENASSN];                     /* |pi.x - pj.x | */

};




/* ================================ */
struct mladen_globals {
/* ================================ */

  FILE *outfilep;
  FILE *oneTimeFlagFilep;
  int called_fluxes;
  long long npart;

  int dump_nr;              /* index of dump */

  struct gizmo_debug_dump * data; /* array for data dump */

} ;




extern struct mladen_globals mladen_globs;


extern void mladen_setup(void);
extern void mladen_cleanup(void);
extern void mladen_dump_after_timestep(void);
extern void mladen_setup_data_dump(long long npart);
extern void mladen_reset_dump_data(void);
extern void mladen_store_particle_data(struct part *p, float h);
// extern void abc(struct part *restrict pi,
//     int pjid,
//     float GSCX,
//     float GSCY,
//     float GSDX,
//     float GSDY,
//     float dwdr,
//     const float r,
//     float hi);

void abc(struct part *restrict pi, 
    int pjid, 
    float GSCX, 
    float GSCY, 
    float GSDX, 
    float GSDY, 
    float dwdr, 
    const float r, 
    float hi){
  /* pi: particle i
   * pjid: id of particle j
   * GSCX, GSCY: Gradient Sum Contribution in x, y direction
   * GSDX, GSDY: Gradient sum dx, dy
   * dwdr: dwdr
   * r: r
   * hi: h_i
   */

  printf("   Storing data for particle %lld\n", pi->id);

  /* first store particle data to be safe */
  mladen_store_particle_data(pi, hi);
  int ind = (int) pi->id;
  struct gizmo_debug_dump * dumploc = &(mladen_globs.data[ind]);
  dumploc->nneigh_grads += 1;
  int ng = dumploc->nneigh_grads;
  if (dumploc->nneigh_grads == MLADENASSN) error("Particle %lld has > %d neighbours\n", pi->id, MLADENASSN);
  dumploc->neighbour_ids_grad[ng] = pjid;
  dumploc->grads_sum_contrib[2*ng] = GSCX;
  dumploc->grads_sum_contrib[2*ng+1] = GSCY;
  dumploc->grads_sum_dx[2*ng] = GSDX;
  dumploc->grads_sum_dx[2*ng+1] = GSDX;
  dumploc->dwdr[ng] = dwdr;
  dumploc->r[ng] = r;
}


#endif /* todo_temporary_globals.h */
