#include "config_maker.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "vector3.h"
#include "utils.h"
#include "math_utils.h"
#include "parameter.h"
#include "system.h"

void createStraightChain(System* system,
	const Parameter* param)
{
	const int32_t num_ptcl = getNumPtcl(param);
	const double blen = getBondLen(param);
	dvec* pos = getPos(system);
	dvec dr = { 0.0, 0.0, 0.0 };
	dr.x = blen;
	clear_dvec(&pos[0]);
	for (int32_t i = 1; i < num_ptcl; i++) {
		pos[i] = add_dvec_new(&pos[i - 1], &dr);
	}
}

void createFlatMesh(System* system,
	const Parameter* param)
{
	const int32_t num_ptcl = getNumPtcl(param);
	if (!isSquareNumber(num_ptcl)) {
		fprintf(stderr, "Error occurs at %s:%d\n", __FILE__, (int32_t)__LINE__);
		fprintf(stderr, "Number of particles should be square number\n");
		fprintf(stderr, "%d is not square number.\n", num_ptcl);
		exit(1);
	}

	dvec* pos = getPos(system);
	const double blen = getBondLen(param);
	const int32_t side_dim = (int32_t)sqrt(num_ptcl);

	int32_t cnt = 0;
	dvec r = { 0.0, 0.0, 0.0 };
	for (int32_t y = 0; y < side_dim; y++) {
		r.y = (y + 0.5) * blen;
		for (int32_t x = 0; x < side_dim; x++) {
			r.x = (x + 0.5) * blen;
			pos[cnt++] = r;
		}
	}
}