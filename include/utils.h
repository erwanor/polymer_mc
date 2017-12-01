#ifndef UTILS_H
#define UTILS_H

#define STR(x) #x
#define CONCATNATE(x, y) x ## y
#define CONCAT(x, y) CONCATNATE(x, y)
#define STRUCT_AT(x, y) CONCAT(x-, >y)
#define UNUSED_PARAMETER(val) (void)(val)
#define UNUSED_FUNCTION(func) (void)(func)

#ifdef DEBUG
#define DEBUG_PRINT(...)                        \
  do {                                          \
    fprintf(stderr, __VA_ARGS__);               \
  } while (0)
#else
#define DEBUG_PRINT(...) do {} while (0)
#endif

#include <stdio.h>
#include <stddef.h>

void* xmalloc(const size_t size);
void xfree(void* ptr);

#include "err_msgs.h"

#define RANGE_CHECK(self, idx)                                          \
  do {                                                                  \
    if ((idx) >= (self)->size) {                                        \
      fprintf(stderr, "%s: at %s %d\n", err_msgs[OUT_OF_RANGE], __FILE__, __LINE__); \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
  } while (0)

#endif
