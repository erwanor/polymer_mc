#ifndef PARAMETER_H
#define PARAMETER_H

#include <stdint.h>
#include <stdbool.h>

#include "vector3.h"
#include "string_c.h"

struct Parameter_t;
typedef struct Parameter_t Parameter;

Parameter* newParameter(const char* path);
void deleteParameter(Parameter* self);

int32_t getNumPtcl(const Parameter* self);
double getBondLen(const Parameter* self);
double getStepLen(const Parameter* self);
int32_t getTotalSteps(const Parameter* self);
int32_t getObserveIntervalMac(const Parameter* self);
int32_t getObserveIntervalMic(const Parameter* self);
uint32_t getRandSeed(const Parameter* self);
double getCfBond(const Parameter* self);
double getCfAngle(const Parameter* self);
const string* getRootDir(const Parameter* self);
const string* getBoundaryName(const Parameter* self);
dvec getBoxlength(const Parameter* self);
int32_t getSideDimx(const Parameter* self);
int32_t getSideDimy(const Parameter* self);

void readParameterFromFile(Parameter* self, const bool is3d);

#endif
