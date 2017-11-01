#ifndef VECTOR3_H
#define VECTOR3_H

#include "utils.h"

typedef struct {
  double x, y, z;
} dvec;

void clear_dvec(dvec* self);
void fill_dvec(dvec* self, const double val);

#define ADD_SUFFIX(a, b) CONCAT(CONCAT(a, _), b)

#define DECL_OP_VEC(OPNAME, VEC) void ADD_SUFFIX(OPNAME, VEC)(VEC* self, const VEC* v0)
#define DECL_OP_VEC_NEW(OPNAME, VEC) VEC CONCAT(ADD_SUFFIX(OPNAME, VEC), _new)(const VEC* v0, const VEC* v1)
#define DECL_OP_SCALAR(OPNAME, VEC, T) void ADD_SUFFIX(OPNAME, scalar)(VEC* self, const T v)
#define DECL_OP_SCALAR_NEW(OPNAME, VEC, T) VEC CONCAT(ADD_SUFFIX(OPNAME, scalar), _new)(const VEC* v0, const T v)

#define DECL_OP_FUNCS(VEC, T, OPNAME, OP) \
DECL_OP_VEC(OPNAME, VEC);\
DECL_OP_VEC_NEW(OPNAME, VEC);\
DECL_OP_SCALAR(OPNAME, VEC, T);\
DECL_OP_SCALAR_NEW(OPNAME, VEC, T);

DECL_OP_FUNCS(dvec, double, add, +)
DECL_OP_FUNCS(dvec, double, sub, -)
DECL_OP_FUNCS(dvec, double, mul, *)
DECL_OP_FUNCS(dvec, double, div, /)

double dvec_dot(const dvec* v0, const dvec* v1);
double norm2(const dvec* self);
double norm(const dvec* self);
void normalize(dvec* vec);

#endif