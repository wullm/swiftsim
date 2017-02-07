#ifndef SWIFT_STRAGGLERS_H
#define SWIFT_STRAGGLERS_H

struct stragglers{
struct spart* stars;
int count;
int size;
};

struct straggler_link{
  struct spart* star;
  struct straggler_link* next;
};

void stragglers_init(struct stragglers* s);
void stragglers_clean(struct stragglers* s);
struct spart* stragglers_add(struct stragglers* s,struct spart* st);

#endif /* SWIFT_STRAGGLERS_H */
