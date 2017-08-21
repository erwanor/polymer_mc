#ifndef TOPOL_H
#define TOPOL_H

#include <stdint.h>
#include "parameter.h"

typedef struct {
	int32_t i0, i1;
} pair;

typedef struct {
	int32_t i0, i1, i2;
} triple;

struct topol_t;
typedef struct topol_t topol;

#if 0
// for 2D mesh
#define NUM_NEIGHBOR_PAIR 4
#define NUM_NEIGHBOR_TRIPLE 12
#else
// for 1D chain
#define NUM_NEIGHBOR_PAIR 2
#define NUM_NEIGHBOR_TRIPLE 3
#endif

typedef struct {
	int32_t num_pair;
	pair pair[NUM_NEIGHBOR_PAIR];
	int32_t num_triple;
	triple triple[NUM_NEIGHBOR_TRIPLE];
} ptclid2topol;

topol* newTopolChain(const Parameter* param);
void deleteTopol(topol* top);

ptclid2topol* newId2Topol(const topol* top, const Parameter* param);
void deleteId2Topol(ptclid2topol* id2top);

int32_t getNumBonds(const topol* top);
const pair* getBondTopol(const topol* top);
int32_t getNumAngles(const topol* top);
const triple* getAngleTopol(const topol* top);

void debugDumpTopolInfo(const topol* top, const Parameter* param);
void debugDumpId2TopolInfo(const ptclid2topol* id2top, const Parameter* param);
#endif