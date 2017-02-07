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
  s->count = 0;
  s->stars = malloc(s->size*sizeof(struct spart));
}

void stragglers_clean(struct stragglers* s){
  
  free(s->stars);
}

struct spart* stragglers_add(struct stragglers* s,struct spart* st){

  if (s->count == s->size)
    error("Not anough space for another star particle!");

  s->stars[s->count] = *st;
  s->count++;
  
  return &s->stars[s->count-1];
}
