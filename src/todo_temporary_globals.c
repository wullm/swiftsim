#include <stdio.h>
#include <string.h>
#include "error.h"
#include "todo_temporary_globals.h"



struct mladen_globals mladen_globs;

/* sets up the necessary stuff. To be called in main.
 * Just drop in everything you need to be done first
 * in here.
 */
void mladen_setup(void){

  mladen_globs.outfilep = fopen("mladen_outputfile_all.txt", "w");  
  mladen_globs.testing = 23;

  message("MESSAGE TO MLADEN: calling init your own global temporary stuff.");
}



void mladen_cleanup(void){
  
  fclose(mladen_globs.outfilep);
}



