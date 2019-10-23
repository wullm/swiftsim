#include <stdio.h>
#include <string.h>
#include "error.h"
#include "todo_temporary_globals.h"


struct mladen_globals mladen_globs;

/* ========================================= */
void mladen_setup(void){
/* ========================================= */
  /* sets up the necessary stuff. To be called in main.
   * Just drop in everything you need to be done first
   * in here.
   */

  mladen_globs.outfilep = fopen("mladen_outputfile_all.txt", "w");  
  mladen_globs.oneTimeFlagFilep = fopen("mladen_flags_all.txt", "w");  
  mladen_globs.called_fluxes = 1;

  message("MESSAGE TO MLADEN: calling init your own global temporary stuff.");


  mladen_globs.dump_nr = 0;

}



/* ========================================= */
void mladen_cleanup(void){
/* ========================================= */
  /* clean up after yourself */
  /* to be called at the end of the run. */
  fclose(mladen_globs.outfilep);
  fclose(mladen_globs.oneTimeFlagFilep);
}







/* ========================================= */
void mladen_setup_data_dump(long long npart){
/* ========================================= */

  /* allocate data dump array, set up necessary stuff.
   * call in main after you established all the necessary
   * parameters necessary. */

  mladen_globs.npart = npart;
  mladen_globs.data = (struct gizmo_debug_dump*) swift_malloc("debug_data", npart*sizeof(struct gizmo_debug_dump));

  /* initialize array to zero */
  mladen_reset_dump_data();

}






/* ======================================= */
void mladen_reset_dump_data(){
/* ======================================= */
  /* reset dump data values */
  struct mladen_globals m = mladen_globs;
  for (int i = 0; i < m.npart; i++){
    struct gizmo_debug_dump *d = &m.data[i];
    d->nneigh = -1;
    d->nneigh_grads = -1;
    d->volume_store = 0;
    d->omega = 0;
    d->wgrads_store[0] = 0;
    d->wgrads_store[1] = 0;
    d->wgrads_store[2] = 0;
    d->pos[0] = 0;
    d->pos[1] = 0;
    d->pos[2] = 0;
    d->h = 0;
    d->id = 0;

    for (int j = 0; j < MLADENASSN; j++){
      d->neighbour_ids[j] = 0;
      d->neighbour_ids_grad[j] = 0;
      d->dwdr[j] = 0;
      d->r[j] = 0;
    }
    for (int j = 0; j < 2*MLADENASSN; j++){
      d->grads_sum_contrib[j] = 0;
      d->Aij[j] = 0;
      d->grads_sum_dx[j] = 0;
    }
  }

  for (int i = 0; i<10; i++){
    printf("Nneigh in data dump is %d\n", mladen_globs.data[i].nneigh);
  }
}







/* ========================================= */
void mladen_dump_after_timestep(void){
/* ========================================= */
  printf("Dumping debugging data %d after timestep\n", mladen_globs.dump_nr);

  /* get filename */
  char filename[200] = "swift-gizmo-debug-dump_";
  char dumpstring[5];
  sprintf(dumpstring, "%04d", mladen_globs.dump_nr);
  strcat(filename, dumpstring);
  strcat(filename, ".txt");
  printf("Filename is %s\n", filename);
  printf("npart is: %lld\n", mladen_globs.npart);

  /* increase dump index */
  mladen_globs.dump_nr += 1;
  /* reset values */
  mladen_reset_dump_data();
}




/* ======================================================== */
void mladen_store_particle_data(struct part *p, float h){
/* ======================================================== */

  /* store default particle data */

  int ind = (int) p->id;
  struct gizmo_debug_dump * dumploc = &(mladen_globs.data[ind]);
  dumploc->id = ind;
  dumploc->pos[0] = p->x[0];
  dumploc->pos[1] = p->x[1];
  dumploc->pos[2] = p->x[2];
  dumploc->h = h;
}


