#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include "vector3.h"
#include "tensor3.h"
#include "boundary.h"

double calcBondEnergy(const dvec* pos0,	const dvec* pos1, const double k, const Boundary* bound);
double calcAngleEnergy(const dvec* pos0, const dvec* pos1, const dvec* pos2, const double k, const Boundary* bound);
dtensor3 calcBondVirial(const dvec* pos0, const dvec* pos1, const double k, const Boundary* bound);
dtensor3 calcAngleVirial(const dvec* pos0, const dvec* pos1, const dvec* pos2, const double k, const Boundary* bound);

#endif