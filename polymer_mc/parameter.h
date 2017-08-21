#ifndef PARAMETER_H
#define PARAMETER_H

#include <stdint.h>
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
unsigned long getRandSeed(const Parameter* self);
double getCfBond(const Parameter* self);
double getCfAngle(const Parameter* self);
string* getRootDir(const Parameter* self);

void readParameterFromFile(Parameter* self);

#endif
