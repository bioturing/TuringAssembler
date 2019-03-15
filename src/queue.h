#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "attribute.h"

typedef struct {
  uint32_t n; //current number of elements
  uint32_t s; //capacity
  uint32_t p; //position of head
  gint_t *e; //array of elements:
} aqueue_t;

aqueue_t *init_aqueue();


void aqueue_add(aqueue_t *q, gint_t e_);

gint_t aqueue_pop(aqueue_t *q);

void aqueue_destroy(aqueue_t *q);
