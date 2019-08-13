#ifndef TODO_TEMPORARY_GLOBALS
#define TODO_TEMPORARY_GLOBALS


struct mladen_globals {

  FILE *outfilep;
  int testing;

} ;

extern struct mladen_globals mladen_globs;


extern void mladen_setup(void);

extern void mladen_cleanup(void);


#endif /* todo_temporary_globals.h */
