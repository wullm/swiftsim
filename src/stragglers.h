#ifndef SWIFT_STRAGGLERS_H
#define SWIFT_STRAGGLERS_H

struct stragglers{
  struct spart* sparts;
  struct gpart* gparts;
  int scount; //Current number of sparts
  int gcount; //Current number of gparts
  int size;
};

struct spart_straggler_link{
  struct spart* sp;
  struct spart_straggler_link* next;
};

struct gpart_straggler_link{
  struct gpart* gp;
  struct gpart_straggler_link* next;
};


void stragglers_init(struct stragglers* s);
void stragglers_clean(struct stragglers* s);
struct spart* stragglers_add_spart(struct stragglers* s, struct spart* sp);
struct gpart* stragglers_add_gpart(struct stragglers* s, struct gpart* gp); 
#endif /* SWIFT_STRAGGLERS_H */
