#ifndef INTERACTIONS_H
#define INTERACTIONS_H

#include "vector3.h"
#include "tensor3.h"
#include "boundary.h"

#include <math.h>

static inline double calcBondEnergy(const dvec* pos0,
                      const dvec* pos1,
                      const double k,
                      const Boundary* bound)
{
  return 0.5 * k * distance2(pos0, pos1, bound);
}

static inline dtensor3 calcBondVirial(const dvec* pos0,
                        const dvec* pos1,
                        const double k,
                        const Boundary* bound)
{
	dvec dr01 = sub_dvec_new(pos1, pos0);
	applyMinimumImageConv(bound, &dr01);
	const dvec dF01 = mul_scalar_new(&dr01, k);
	return dtensor3_dot(&dF01, &dr01);
}

static inline double calcAngleEnergy(const dvec* pos0,
                       const dvec* pos1,
                       const dvec* pos2,
                       const double k,
                       const Boundary* bound)
{
	return k * (1.0 - cos_angle(pos0, pos1, pos2, bound));
}

static inline dtensor3 calcAngleVirial(const dvec* pos0,
	const dvec* pos1,
	const dvec* pos2,
	const double k,
	const Boundary* bound)
{
	dvec dr10 = sub_dvec_new(pos0, pos1); // 1 -> 0
	dvec dr12 = sub_dvec_new(pos2, pos1); // 1 -> 2
	applyMinimumImageConv(bound, &dr10);
	applyMinimumImageConv(bound, &dr12);

	const double dr10_norm2 = norm2(&dr10);
	const double dr12_norm2 = norm2(&dr12);

	const double dr10_norm = sqrt(dr10_norm2);
	const double dr12_norm = sqrt(dr12_norm2);

	double cs = dvec_dot(&dr10, &dr12) / (dr10_norm * dr12_norm);
	if (cs > 1.0) cs = 1.0;
	if (cs < -1.0) cs = -1.0;

	const double a11 = k * cs / dr10_norm2;
	const double a12 = -k / (dr10_norm * dr12_norm);
	const double a22 = k * cs / dr12_norm2;

	const dvec dF0_0 = mul_scalar_new(&dr10, a11);
	const dvec dF0_1 = mul_scalar_new(&dr12, a12);
	const dvec dF0 = add_dvec_new(&dF0_0, &dF0_1);
	const dtensor3 dF0_dr10 = dtensor3_dot(&dF0, &dr10);

	const dvec dF1_0 = mul_scalar_new(&dr12, a22);
	const dvec dF1_1 = mul_scalar_new(&dr10, a12);
	const dvec dF1 = add_dvec_new(&dF1_0, &dF1_1);
	const dtensor3 dF1_dr12 = dtensor3_dot(&dF1, &dr12);

	return dtensor3_add_new(&dF0_dr10, &dF1_dr12);
}

#endif
