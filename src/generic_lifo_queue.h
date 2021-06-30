//
// Created by yuyttenh on 30/06/2021.
//

/* QUEUE_NAME and QUEUE_TYPE must be defined before importing this header! */

#define PASTE(x, y) x##_##y

#ifndef QUEUE_NAME
#define _NAME(f) PASTE(f, lifo_queue)
#define QUEUE_NAME _NAME(QUEUE_TYPE)
#endif

#define _INIT(f) PASTE(f, init)
#define QUEUE_INIT _INIT(QUEUE_NAME)

#define _DESTROY(f) PASTE(f, destroy)
#define QUEUE_DESTROY _DESTROY(QUEUE_NAME)

#define _RESET(f) PASTE(f, reset)
#define QUEUE_RESET _RESET(QUEUE_NAME)

#define _IS_EMPTY(f) PASTE(f, is_empty)
#define QUEUE_IS_EMPTY _IS_EMPTY(QUEUE_NAME)

#define _PUSH(f) PASTE(f, push)
#define QUEUE_PUSH _PUSH(QUEUE_NAME)

#define _POP(f) PASTE(f, pop)
#define QUEUE_POP _POP(QUEUE_NAME)

struct QUEUE_NAME {
  QUEUE_TYPE *values;
  int size;
  int index;
};

inline static void QUEUE_INIT(struct QUEUE_NAME *q, int size) {
  q->values = (QUEUE_TYPE *)malloc(size * sizeof(QUEUE_TYPE));
  q->size = size;
  q->index = 0;
}

inline static void QUEUE_DESTROY(struct QUEUE_NAME *q) {
  free(q->values);
}

inline static void QUEUE_RESET(struct QUEUE_NAME *q) {
  q->index = 0;
}

inline static int QUEUE_IS_EMPTY(struct QUEUE_NAME *q) {
  return q->index == 0;
}

inline static void QUEUE_PUSH(struct QUEUE_NAME *q, QUEUE_TYPE value) {
  if (q->size == q->index) {
    q->size <<= 1;
    q->values = realloc(q->values, q->size * sizeof(QUEUE_TYPE));
  }
  q->values[q->index] = value;
  q->index++;
}

inline static QUEUE_TYPE QUEUE_POP(struct QUEUE_NAME *q) {
#ifdef QUEUE_SAFETY_CHECKS
  if (QUEUE_IS_EMPTY(q)) {
    fprintf(stderr, "Trying to pop from empty queue!");
    abort();
  }
#endif
  q->index--;
  return q->values[q->index];
}

/* Undefine QUEUE_NAME and QUEUE_TYPE */
#undef QUEUE_NAME
#undef QUEUE_TYPE