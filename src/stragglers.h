#ifndef SWIFT_STRAGGLERS_H
#define SWIFT_STRAGGLERS_H

struct stragglers{
  struct part* parts;
  struct gpart* gparts;
  struct spart* sparts;
  struct xpart* xparts;
  int count; //Current number of parts
  int scount; //Current number of sparts
  int gcount; //Current number of gparts
  int xcount; //Current number of xparts, should always be the same as count
  int size;
};

struct straggler_link{
  struct part* p;
  struct straggler_link* next;
};

struct g_straggler_link{
  struct gpart* gp;
  struct g_straggler_link* next;
};

struct s_straggler_link{
  struct spart* sp;
  struct s_straggler_link* next;
};

struct x_straggler_link{
  struct xpart* xp;
  struct x_straggler_link* next;
};

void stragglers_init(struct stragglers* s);
void stragglers_clean(struct stragglers* s);

struct part* stragglers_add_part(struct stragglers* s, struct part* p);
struct gpart* stragglers_add_gpart(struct stragglers* s, struct gpart* gp); 
struct spart* stragglers_add_spart(struct stragglers* s, struct spart* sp);
struct xpart* stragglers_add_xpart(struct stragglers* s, struct xpart* xp);
#endif /* SWIFT_STRAGGLERS_H */
