#ifndef SYSTEM_H
#define SYSTEM_H

#include "vector3.h"
#include "parameter.h"
#include "topol.h"

struct System_t;
typedef struct System_t System;

typedef void(*confMaker)(System* system, const Parameter* param);
typedef topol*(*topolMaker)(const Parameter*);

System* newSystem(void);
void deleteSystem(System* self);

topol* getTopol(const System* self);
ptclid2topol* getPtclId2Topol(const System* self);
dvec* getPos(const System* self);
double getAcceptRatio(const System* self);

void initializeSystem(System* self, const Parameter* param,	confMaker conf_make, topolMaker topol_make);
void executeSimulation(System* self, const Parameter* param);

void readRestartConfig(System* self, const Parameter* param);
void writeFinalConfig(System* self, const Parameter* param);

#endif