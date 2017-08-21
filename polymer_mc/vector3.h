#ifndef VECTOR3_H
#define VECTOR3_H

typedef struct {
  double x, y, z;
} dvec;

void clear_dvec(dvec* self);
void fill_dvec(dvec* self, const double val);

double get_x(const dvec* self);
void set_x(dvec* self, const double x);
double get_y(const dvec* self);
void set_y(dvec* self, const double y);
double get_z(const dvec* self);
void set_z(dvec* self, const double z);

#define CONCATENATE(a, b) a ## b
#define CONCAT(a, b) CONCATENATE(a, b)
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
DECL_OP_FUNCS(dvec, double, div, / )

double dvec_dot(const dvec* v0, const dvec* v1);
double norm2(const dvec* self);
double norm(const dvec* self);
double distance2(const dvec* pos0, const dvec* pos1);
double distance(const dvec* pos0, const dvec* pos1);
double cos_angle(const dvec* pos0, const dvec* pos1, const dvec* pos2);
void normalize(dvec* vec);

#endif