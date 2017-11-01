#ifndef SYSTEM_H
#define SYSTEM_H

#include "vector3.h"

struct topol_t;
typedef struct topol_t topol;

struct ptclid2topol_t;
typedef struct ptclid2topol_t ptclid2topol;

struct System_t;
typedef struct System_t System;

struct Parameter_t;
typedef struct Parameter_t Parameter;

struct Boundary_t;
typedef struct Boundary_t Boundary;

typedef void(*confMaker)(System* system, const Parameter* param);
typedef topol*(*topolMaker)(const Parameter*, const Boundary* bound);

System* newSystem(void);
void deleteSystem(System* self);

topol* getTopol(const System* self);
ptclid2topol* getPtclId2Topol(const System* self);
dvec* getPos(const System* self);
double getAcceptRatio(const System* self);

void initializeSystem(System* self, const Boundary* boundary, const Parameter* param, confMaker conf_make, topolMaker topol_make);
void executeSimulation(System* self, const Boundary* boundary, const Parameter* param);

void readRestartConfig(System* self, const Parameter* param);
void writeFinalConfig(System* self, const Parameter* param);

#endif