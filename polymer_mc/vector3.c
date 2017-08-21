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

double get_x(const dvec* self)
{
	return self->x;
}

void set_x(dvec* self, const double x)
{
	self->x = x;
}

double get_y(const dvec* self)
{
	return self->y;
}

void set_y(dvec* self, const double y)
{
	self->y = y;
}

double get_z(const dvec* self)
{
	return self->z;
}

void set_z(dvec* self, const double z)
{
	self->z = z;
}

#define DEFINE_OP_FUNCS(VEC, T, OPNAME, OP) \
DECL_OP_VEC(OPNAME, VEC) \
{\
	self->x CONCAT(OP, =) v0->x;\
	self->y CONCAT(OP, =) v0->y;\
	self->z CONCAT(OP, =) v0->z;\
}\
\
DECL_OP_VEC_NEW(OPNAME, VEC) \
{\
	VEC ret; \
	ret.x = v0->x OP v1->x;\
	ret.y = v0->y OP v1->y;\
	ret.z = v0->z OP v1->z;\
	return ret;\
}\
\
DECL_OP_SCALAR(OPNAME, VEC, T)\
{\
	self->x CONCAT(OP, =) v; \
	self->y CONCAT(OP, =) v; \
	self->z CONCAT(OP, =) v; \
}\
\
DECL_OP_SCALAR_NEW(OPNAME, VEC, T)\
{\
	VEC ret;\
	ret.x = v0->x OP v;\
	ret.y = v0->y OP v;\
	ret.z = v0->z OP v;\
	return ret; \
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

// NOTE: return dr2
double norm2(const dvec* self)
{
	return dvec_dot(self, self);
}

// NOTE: return sqrt(dr2)
double norm(const dvec* self)
{
	return sqrt(norm2(self));
}

double distance2(const dvec* pos0,
	const dvec* pos1)
{
	const dvec dr = sub_dvec_new(pos0, pos1);
	return dvec_dot(&dr, &dr);
}

double distance(const dvec* pos0,
	const dvec* pos1)
{
	return sqrt(distance2(pos0, pos1));
}

// NOTE: return (dr01*dr21) / (|dr01|*|dr21|)
double cos_angle(const dvec* pos0,
	const dvec* pos1,
	const dvec* pos2)
{
	const dvec dr01 = sub_dvec_new(pos1, pos0); // 0 -> 1
	const dvec dr12 = sub_dvec_new(pos2, pos1); // 1 -> 2
	const double dr01_dr12 = dvec_dot(&dr01, &dr12);
	const double dr01_norm = norm(&dr01);
	const double dr12_norm = norm(&dr12);
	return dr01_dr12 / (dr01_norm * dr12_norm);
}

void normalize(dvec* vec)
{
	div_scalar(vec, norm2(vec));
}

void rotateAlongZ(dvec* vec, const double theta)
{
	const double new_vec_x = vec->x * cos(theta) - vec->y * sin(theta);
	const double new_vec_y = vec->x * sin(theta) + vec->y * cos(theta);
	vec->x = new_vec_x;
	vec->y = new_vec_y;
}

/*void rotateAlongAxis(dvec* vec, const dvec* axis)
{

}*/