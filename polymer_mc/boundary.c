#include "boundary.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "utils.h"
#include "system.h"

struct Boundary_t {
	dvec box_leng;
	dvec hbox_leng;
	BOUNDARY_TYPE type;
};

// NOTE: boundary origin is {0.0, 0.0, 0.0}
Boundary* newBoundary(const string* type_name)
{
	Boundary* bound = (Boundary*)xmalloc(sizeof(Boundary));
	bound->type = getBoundaryTypeFromName(type_name);
	return bound;
}

void deleteBoundary(Boundary* bound)
{
	xfree(bound);
}

BOUNDARY_TYPE getBoundaryType(const Boundary* bound)
{
	return bound->type;
}

const char* getBoundaryNameFromType(BOUNDARY_TYPE type)
{
	switch (type) {
	case FREE:
		return "free";
	case PERIODIC:
		return "periodic";
	default:
		fprintf(stderr, "Unknown type\n");
		return NULL;
	}
}

#define COMPARE_BOUNDARY_TYPE(boundary_name, BTYPE) \
(0 == strncmp(string_to_char(boundary_name), getBoundaryNameFromType(BTYPE), strlen(getBoundaryNameFromType(BTYPE))))

BOUNDARY_TYPE getBoundaryTypeFromName(const string* boundary_name)
{
	if (COMPARE_BOUNDARY_TYPE(boundary_name, FREE)) {
		return FREE;
	} else if (COMPARE_BOUNDARY_TYPE(boundary_name, PERIODIC)) {
		return PERIODIC;
	} else {
		fprintf(stderr, "Error occurs in getBoundaryTypeFromName at %s:%d\n", __FILE__, (int32_t)__LINE__);
		fprintf(stderr, "Unknown boundary name %s\n", string_to_char(boundary_name));
		exit(1);
	}
}

void setBoxLength(Boundary* self,
	const dvec length)
{
	if (self->type == FREE) {
		fprintf(stderr, "cannot set box length for %s\n", getBoundaryNameFromType(self->type));
		exit(1);
	}
	self->box_leng  = length;
	self->hbox_leng = mul_scalar_new(&length, 0.5);
}

void applyBoundaryCond(const Boundary* self,
	dvec* pos)
{
	if (self->type == FREE) return;
	if (pos->x < 0.0) pos->x += self->box_leng.x;
	if (pos->x >= self->box_leng.x) pos->x -= self->box_leng.x;
	if (pos->y < 0.0) pos->y += self->box_leng.y;
	if (pos->y >= self->box_leng.y) pos->y -= self->box_leng.y;
	if (pos->z < 0.0) pos->z += self->box_leng.z;
	if (pos->z >= self->box_leng.z) pos->z -= self->box_leng.z;
}

void applyBoundaryCondForSystem(const Boundary* self,
	System* system,
	const Parameter* param)
{
	if (self->type == FREE) return;
	const int32_t num_ptcls = getNumPtcl(param);
	dvec* pos = getPos(system);
	for (int32_t i = 0; i < num_ptcls; i++) {
		applyBoundaryCond(self, &pos[i]);
	}
}

void applyMinimumImageConv(const Boundary* self,
	dvec* dr)
{
	if (self->type == FREE) return;
	if (dr->x >  self->hbox_leng.x) dr->x -= self->box_leng.x;
	if (dr->x < -self->hbox_leng.x) dr->x += self->box_leng.x;
	if (dr->y >  self->hbox_leng.y) dr->y -= self->box_leng.y;
	if (dr->y < -self->hbox_leng.y) dr->y += self->box_leng.y;
	if (dr->z >  self->hbox_leng.z) dr->z -= self->box_leng.z;
	if (dr->z < -self->hbox_leng.z) dr->z += self->box_leng.z;
}

double distance2(const dvec* pos0,
	const dvec* pos1,
	const Boundary* bound)
{
	dvec dr = sub_dvec_new(pos0, pos1);
	applyMinimumImageConv(bound, &dr);
	return dvec_dot(&dr, &dr);
}

double distance(const dvec* pos0,
	const dvec* pos1,
	const Boundary* bound)
{
	return sqrt(distance2(pos0, pos1, bound));
}

// NOTE: return (dr01*dr21) / (|dr01|*|dr21|)
double cos_angle(const dvec* pos0,
	const dvec* pos1,
	const dvec* pos2,
	const Boundary* bound)
{
	dvec dr01 = sub_dvec_new(pos1, pos0); // 0 -> 1
	dvec dr12 = sub_dvec_new(pos2, pos1); // 1 -> 2
	applyMinimumImageConv(bound, &dr01);
	applyMinimumImageConv(bound, &dr12);
	const double dr01_dr12 = dvec_dot(&dr01, &dr12);
	const double dr01_norm = norm(&dr01);
	const double dr12_norm = norm(&dr12);
	return dr01_dr12 / (dr01_norm * dr12_norm);
}