#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include "attribute.h"
#include "verbose.h"
#include "queue.h"

#define INIT_SIZE 100

aqueue_t *init_aqueue()
{
  aqueue_t *q = (aqueue_t *)malloc(sizeof(aqueue_t));
  if (q == NULL)
    __ERROR("Cannot allocate memory\n");
  
  q->e = (gint_t *)calloc(INIT_SIZE, sizeof(gint_t));
  if (q->e == NULL)
    __ERROR("Cannot allocate memory\n");
  q->s = INIT_SIZE;
  q->p = 0;
  q->n = 0;
  return q;
}

void aqueue_add(aqueue_t *q, gint_t e_)
{
   if (q->p + q->n == q->s){
    if (q->p > 0)
      memmove(q->e, q->e + q->p, q->n * sizeof(gint_t));
    else{
      q->e = (gint_t *)realloc(q->e, sizeof(gint_t) * (q->s << 1));
      q->s <<= 1;
    }
    q->p = 0;
   } 
   q->e[q->p + q->n] = e_;
   q->n++;
}

gint_t aqueue_pop(aqueue_t *q)
{
  assert(q->n > 0);
  gint_t ret;
  ret = q->e[q->p];
  if (q->n == 1){
    q->p = 0;
    q->n = 0;
  }
  else {
    q->p++;
    q->n--;
  }
  return ret;
}

void aqueue_destroy(aqueue_t *q)
{
  assert(q != NULL);
  assert(q->e != NULL);
  free(q->e);
  free(q);
}

