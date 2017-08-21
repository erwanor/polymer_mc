#include "config_maker.h"

#include <stdio.h>

#include "vector3.h"

void createStraightChain(System* system,
	const Parameter* param)
{
	const int32_t num_ptcl = getNumPtcl(param);
	const double blen = getBondLen(param);
	dvec* pos = getPos(system);
	dvec dr = {0.0};
	dr.x = blen;
	clear_dvec(&pos[0]);
	for (int32_t i = 1; i < num_ptcl; i++) {
		pos[i] = add_dvec_new(&pos[i - 1], &dr);
	}
}

/*void createFlatMesh(System* system,
	const Parameter* param)
{

}*/