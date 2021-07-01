//
// Created by yuyttenh on 01/07/2021.
//

#ifndef CVORONOI_TEST_QUEUES_C
#define CVORONOI_TEST_QUEUES_C

#include <assert.h>
#include <stdio.h>

#include "queues.h"
#include "tuples.h"

inline static void test_int_lifo_queue(){
  struct int_lifo_queue q;
  int_lifo_queue_init(&q, 1);
  assert(int_lifo_queue_is_empty(&q));

  int_lifo_queue_push(&q, 2);
  int_lifo_queue_push(&q, 3);
  int_lifo_queue_push(&q, 31);
  assert(q.index == 3);
  assert(q.size == 4);

  int v = int_lifo_queue_pop(&q);
  assert(v == 31);
  assert(q.index == 2);
  v = int_lifo_queue_pop(&q);
  assert(v == 3);
  assert(q.index == 1);
  v = int_lifo_queue_pop(&q);
  assert(v == 2);
  assert(int_lifo_queue_is_empty(&q));

  int_lifo_queue_destroy(&q);
}

inline static void test_int3_fifo_queue() {
  struct int3_fifo_queue q;
  int3_fifo_queue_init(&q, 1);
  assert(int3_fifo_queue_is_empty(&q));

  int3_fifo_queue_push(&q, (int3){._0 = 2, ._1 = 1, ._2 = 3});
  int3_fifo_queue_push(&q, (int3){._0 = 3, ._1 = 2, ._2 = 1});
  int3_fifo_queue_push(&q, (int3){._0 = 21, ._1 = 11, ._2 = 31});
  assert(q.end == 3);
  assert(q.start == 0);
  assert(q.size == 4);

  int3 v = int3_fifo_queue_pop(&q);
  assert(v._0 == 2 && v._1 == 1 && v._2 == 3);
  assert(q.start == 1);
  v = int3_fifo_queue_pop(&q);
  assert(v._0 == 3 && v._1 == 2 && v._2 == 1);
  assert(q.start == 2);
  v = int3_fifo_queue_pop(&q);
  assert(v._0 == 21 && v._1 == 11 && v._2 == 31);
  assert(int3_fifo_queue_is_empty(&q));

  int3_fifo_queue_destroy(&q);
}

int main() {
  test_int_lifo_queue();

  test_int3_fifo_queue();

  printf("Succes!");
}

#endif  // CVORONOI_TEST_QUEUES_C
