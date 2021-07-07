//
// Created by yuyttenh on 30/06/2021.
//

/**
 * @file queues.h
 *
 * @brief Generates code for a int LIFO queue and an int3 FIFO queue
 */

#ifndef CVORONOI_QUEUES_H
#define CVORONOI_QUEUES_H

#include "tuples.h"

#define QUEUE_SAFETY_CHECKS

#define QUEUE_TYPE int
#include "generic_lifo_queue.h"

#define QUEUE_TYPE int3
#include "generic_fifo_queue.h"

#endif  // CVORONOI_QUEUES_H
