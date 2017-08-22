#ifndef OBSERBER_H
#define OBSERBER_H

#include <stdint.h>

#include "string_c.h"
#include "system.h"
#include "parameter.h"
#include "boundary.h"

struct Observer_t;
typedef struct Observer_t Observer;

Observer* newObserver(const string* dir_path);
void deleteObserver(Observer* self);

void observeMicroVars(Observer* self, const int32_t mc_steps, const System* system, const Boundary* bound, const Parameter* param);
void observeMacroVars(Observer* self, const int32_t mc_steps, const System* system, const Boundary* bound, const Parameter* param);

#endif