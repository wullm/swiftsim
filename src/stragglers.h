#ifndef SWIFT_STRAGGLERS_H
#define SWIFT_STRAGGLERS_H

struct stragglers{
  struct gpart* gparts;
  struct spart* sparts;
  int scount; //Current number of sparts
  int gcount; //Current number of gparts
  int size;
};

struct g_straggler_link{
  struct gpart* gp;
  struct g_straggler_link* next;
};

struct s_straggler_link{
  struct spart* sp;
  struct s_straggler_link* next;
};

void stragglers_init(struct stragglers* s);
void stragglers_clean(struct stragglers* s);

struct gpart* stragglers_add_gpart(struct stragglers* s, struct gpart* gp); 
struct spart* stragglers_add_spart(struct stragglers* s, struct spart* sp);
#endif /* SWIFT_STRAGGLERS_H */
