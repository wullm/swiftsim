/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "stragglers.h"

/* Local headers. */
#include "atomic.h"
#include "const.h"
#include "cooling.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "lock.h"
#include "memswap.h"
#include "minmax.h"
#include "runner.h"
#include "stars.h"
#include "threadpool.h"
#include "tools.h"
void stragglers_init(struct stragglers* s){

  s->size = 1000;
  s->scount = 0;
  s->gcount = 0;
  s->gparts = malloc(s->size*sizeof(struct gpart));
  s->sparts = malloc(s->size*sizeof(struct spart));  
}

void stragglers_clean(struct stragglers* s){

  free(s->gparts);
  free(s->sparts);
}

struct gpart* stragglers_add_gpart(struct stragglers* s,struct gpart* gp){

  if (s->gcount == s->size)
    error("Not anough space for another straggler gpart!");

  s->gparts[s->gcount] = *gp;
  s->gcount++;
  
  return &s->gparts[s->gcount-1];
}

struct spart* stragglers_add_spart(struct stragglers* s,struct spart* sp){

  if (s->scount == s->size)
    error("Not anough space for another straggler spart!");

  s->sparts[s->scount] = *sp;
  s->scount++;
  
  return &s->sparts[s->scount-1];
}

