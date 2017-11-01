#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "vector3.h"
#include "parameter.h"

struct System_t;
typedef struct System_t System;

typedef enum {
  PERIODIC = 0,
  FREE,
} BOUNDARY_TYPE;

struct Boundary_t;
typedef struct Boundary_t Boundary;

Boundary* newBoundary(const string* boundary_name);
void deleteBoundary(Boundary* bound);

void setBoxLength(Boundary* self, const dvec length);
void applyBoundaryCond(const Boundary* self, dvec* pos);
void applyBoundaryCondForSystem(const Boundary* self, System* system, const Parameter* param);
void applyMinimumImageConv(const Boundary* self, dvec* dr);

BOUNDARY_TYPE getBoundaryType(const Boundary* bound);
const char* getBoundaryNameFromType(BOUNDARY_TYPE type);
BOUNDARY_TYPE getBoundaryTypeFromName(const string* boundary_name);
double distance2(const dvec* pos0, const dvec* pos1, const Boundary* bound);
double distance(const dvec* pos0, const dvec* pos1, const Boundary* bound);
double cos_angle(const dvec* pos0, const dvec* pos1, const dvec* pos2, const Boundary* bound);

#endif