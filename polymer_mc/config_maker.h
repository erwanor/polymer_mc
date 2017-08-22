#ifndef CONFIG_MAKER_H
#define CONFIG_MAKER_H

struct System_t;
typedef struct System_t System;

struct Parameter_t;
typedef struct Parameter_t Parameter;

void createStraightChain(System* system, const Parameter* param);
void createFlatMesh(System* system, const Parameter* param);

#endif