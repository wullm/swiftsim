#include "active.h"


/**
 * @brief Calculate the volume interaction between particle i and particle j
 *
 * The volume is in essence the same as the weighted number of neighbours in a
 * classical SPH density calculation.
 *
 * We also calculate the components of the matrix E, which is used for second
 * order accurate gradient calculations and for the calculation of the interface
 * surface areas.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
void runner_iact_density(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

  /* Get r and h inverse. */
  const float r = sqrtf(r2);

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);
#ifdef WITH_IVANOVA
  /* TODO: don't forget about the nonsym part! */
  const float hidp1 = pow_dimension_plus_one(hi_inv);
  pi->density.wgrads[0] += hidp1 * wi_dx * dx[0] / r;
  pi->density.wgrads[1] += hidp1 * wi_dx * dx[1] / r;
  pi->density.wgrads[2] += hidp1 * wi_dx * dx[2] / r;

  // TODO: temporary
  mladen_store_neighbour_data(
        /* pi    */ pi,
        /* pj ID */ pj->id,
        /* GSCX= */ hidp1 * wi_dx * dx[0]/r,
        /* GSCY= */ hidp1 * wi_dx * dx[1]/r,
        /* GSDX= */ dx[0],
        /* GSDY= */ dx[1],
        /* dwdr= */ hidp1 * wi_dx,
        /* r=    */ r,
        /* hi =  */ hi );
#endif

  /* these are eqns. (1) and (2) in the summary */
  pi->geometry.volume += wi;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

  pi->geometry.centroid[0] -= dx[0] * wi;
  pi->geometry.centroid[1] -= dx[1] * wi;
  pi->geometry.centroid[2] -= dx[2] * wi;

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + xj * wj_dx);
#ifdef WITH_IVANOVA
  const float hjdp1 = pow_dimension_plus_one(hj_inv);
  pj->density.wgrads[0] -= hjdp1 * wj_dx * dx[0] / r;
  pj->density.wgrads[1] -= hjdp1 * wj_dx * dx[1] / r;
  pj->density.wgrads[2] -= hjdp1 * wj_dx * dx[2] / r;

  // TODO: temporary
  mladen_store_neighbour_data(pj, pi->id,
        /* GSCX= */ -hjdp1 * wj_dx * dx[0]/r,
        /* GSCY= */ -hjdp1 * wj_dx * dx[1]/r,
        /* GSDX= */ -dx[0],
        /* GSDY= */ -dx[1],
        /* dwdr= */ hjdp1 * wj_dx,
        /* r= */ r, hj );
#endif

  /* these are eqns. (1) and (2) in the summary */
  pj->geometry.volume += wj;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pj->geometry.matrix_E[k][l] += dx[k] * dx[l] * wj;

  pj->geometry.centroid[0] += dx[0] * wj;
  pj->geometry.centroid[1] += dx[1] * wj;
  pj->geometry.centroid[2] += dx[2] * wj;
}





/**
 * @brief Calculate the volume interaction between particle i and particle j:
 * non-symmetric version
 *
 * The volume is in essence the same as the weighted number of neighbours in a
 * classical SPH density calculation.
 *
 * We also calculate the components of the matrix E, which is used for second
 * order accurate gradient calculations and for the calculation of the interface
 * surface areas.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
void runner_iact_nonsym_density(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

  float wi, wi_dx;

  /* Get r and h inverse. */
  const float r = sqrtf(r2);

  const float hi_inv = 1.f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);
#ifdef WITH_IVANOVA
  /* TODO: don't forget about the sym part! */
  const float hidp1 = pow_dimension_plus_one(hi_inv);
  pi->density.wgrads[0] += hidp1 * wi_dx * dx[0] / r;
  pi->density.wgrads[1] += hidp1 * wi_dx * dx[1] / r;
  pi->density.wgrads[2] += hidp1 * wi_dx * dx[2] / r;
    
  // TODO: temporary
  mladen_store_neighbour_data(
        /* pi    */ pi,
        /* pj ID */ pj->id,
        /* GSCX= */ hidp1 * wi_dx * dx[0]/r,
        /* GSCY= */ hidp1 * wi_dx * dx[1]/r,
        /* GSDX= */ dx[0],
        /* GSDY= */ dx[1],
        /* dwdr= */ hidp1 * wi_dx,
        /* r=    */ r,
        /* hi =  */ hi );
#endif

  pi->geometry.volume += wi;
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      pi->geometry.matrix_E[k][l] += dx[k] * dx[l] * wi;

  pi->geometry.centroid[0] -= dx[0] * wi;
  pi->geometry.centroid[1] -= dx[1] * wi;
  pi->geometry.centroid[2] -= dx[2] * wi;
}



