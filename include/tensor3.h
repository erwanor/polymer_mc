#ifndef TENSOR3_H
#define TENSOR3_H

#include "vector3.h"

typedef struct dtensor3_t {
  double xx, xy, xz;
  double yx, yy, yz;
  double zx, zy, zz;
} dtensor3;

dtensor3 dtensor3_dot(const dvec* v0, const dvec* v1);
void dtensor3_add(dtensor3* t0, const dtensor3* t1);
dtensor3 dtensor3_add_new(const dtensor3* t0, const dtensor3* t1);
void dtensor3_sub(dtensor3* t0, const dtensor3* t1);
dtensor3 dtensor3_sub_new(const dtensor3* t0, const dtensor3* t1);
void dtensor3_clear(dtensor3* t);
void dtensor3_mul_scalar(dtensor3* t, const double k);

#endif
