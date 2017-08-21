#ifndef EVOLVER_H
#define EVOLVER_H

#include <stdbool.h>

#include "system.h"
#include "parameter.h"
#include "mt_rand.h"
#include "vector3.h"

dvec kickParticle(const dvec* pos, const double disp, MTstate* mtst);
double calcBondEnergy(const dvec* pos0, const dvec* pos1, const double k);
double calcAngleEnergy(const dvec* pos0, const dvec* pos1, const dvec* pos2, const double k);
bool newStateIsAccepted(const double deltaE, MTstate* mtst);
double evolveMc(System* system, const Parameter* param, MTstate* mtst);

#endif