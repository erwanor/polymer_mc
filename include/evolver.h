#ifndef EVOLVER_H
#define EVOLVER_H

#include <stdbool.h>
#include "vector3.h"
#include "boundary.h"

struct System_t;
typedef struct System_t System;

struct Parameter_t;
typedef struct Parameter_t Parameter;

struct Boundary_t;
typedef struct Boundary_t Boundary;

struct MTstate_t;
typedef struct MTstate_t MTstate;

double evolveMc(System *system, const Parameter *param, const Boundary *bound, MTstate *mtst);
double boundaryDistance(const dvec pos1, const dvec pos2, const Boundary *bound);
bool particuleBoundary(dvec point, const Boundary *bound);
bool checkParticleOverlap(const dvec *os, int32_t num_ptcl, const Boundary *bound);
#endif