#include "utils.h"

#include <stdlib.h>

void* xmalloc(const size_t size) {
  void* ptr_ret = malloc(size);
  if (!ptr_ret) {
    perror("xmalloc");
    exit(EXIT_FAILURE);
  }
  return ptr_ret;
}

void xfree(void* ptr) {
  free(ptr);
}
