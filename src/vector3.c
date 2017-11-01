#include "vector3.h"

#include <math.h>

void clear_dvec(dvec* self)
{
  fill_dvec(self, 0.0);
}

void fill_dvec(dvec* self, const double val)
{
  self->x = self->y = self->z = val;
}

#define DEFINE_OP_FUNCS(VEC, T, OPNAME, OP)     \
  DECL_OP_VEC(OPNAME, VEC)                      \
  {                                             \
    self->x CONCAT(OP, =) v0->x;                \
    self->y CONCAT(OP, =) v0->y;                \
    self->z CONCAT(OP, =) v0->z;                \
  }                                             \
                                                \
  DECL_OP_VEC_NEW(OPNAME, VEC)                  \
  {                                             \
    VEC ret;                                    \
    ret.x = v0->x OP v1->x;                     \
    ret.y = v0->y OP v1->y;                     \
    ret.z = v0->z OP v1->z;                     \
    return ret;                                 \
  }                                             \
                                                \
  DECL_OP_SCALAR(OPNAME, VEC, T)                \
  {                                             \
    self->x CONCAT(OP, =) v;                    \
    self->y CONCAT(OP, =) v;                    \
    self->z CONCAT(OP, =) v;                    \
  }                                             \
                                                \
  DECL_OP_SCALAR_NEW(OPNAME, VEC, T)            \
  {                                             \
    VEC ret;                                    \
    ret.x = v0->x OP v;                         \
    ret.y = v0->y OP v;                         \
    ret.z = v0->z OP v;                         \
    return ret;                                 \
  }

DEFINE_OP_FUNCS(dvec, double, add, +)
DEFINE_OP_FUNCS(dvec, double, sub, -)
DEFINE_OP_FUNCS(dvec, double, mul, *)
DEFINE_OP_FUNCS(dvec, double, div, /)

double dvec_dot(const dvec* v0,
                const dvec* v1)
{
  return v0->x * v1->x + v0->y * v1->y + v0->z * v1->z;
}

double norm2(const dvec* self)
{
  return dvec_dot(self, self);
}

double norm(const dvec* self)
{
  return sqrt(norm2(self));
}

void normalize(dvec* vec)
{
  div_scalar(vec, norm2(vec));
}
