#ifndef EVOLVER_H
#define EVOLVER_H

#include <stdbool.h>

struct System_t;
typedef struct System_t System;

struct Parameter_t;
typedef struct Parameter_t Parameter;

struct Boundary_t;
typedef struct Boundary_t Boundary;

struct MTstate_t;
typedef struct MTstate_t MTstate;

double evolveMc(System* system, const Parameter* param, const Boundary* bound, MTstate* mtst);

#endif