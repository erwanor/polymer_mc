#include "parameter.h"

#include <assert.h>
#include <math.h>

#include "utils.h"

struct Parameter_t
{
	int32_t num_ptcl;
	int32_t total_steps;
	int32_t observe_interval_mic;
	int32_t observe_interval_mac;
	string* root_dir;
	double bond_len;
	double step_len;
	double cf_bond;
	double cf_angle;
	unsigned long rand_seed;
};

static void initializeParameter(Parameter* self);
static void checkParameterIsLoaded(const Parameter* self);

Parameter* newParameter(const char* path)
{
	Parameter* self = (Parameter*)xmalloc(sizeof(Parameter));
	initializeParameter(self);
	self->root_dir = new_string_from_char(path);
	return self;
}

void deleteParameter(Parameter* self)
{
	delete_string(self->root_dir);
	xfree(self);
}

static void initializeParameter(Parameter* self)
{
	self->num_ptcl = -1;
	self->bond_len = -1;
	self->step_len = -1;
	self->cf_bond = nan("");
	self->cf_angle = nan("");
	self->total_steps = -1;
	self->observe_interval_mac = -1;
	self->observe_interval_mic = -1;
	self->rand_seed = 0xffffffff;
	self->root_dir = NULL;
}

static void checkParameterIsLoaded(const Parameter* self)
{
	assert(self->num_ptcl > 0);
	assert(self->bond_len > 0);
	assert(self->step_len > 0);
	assert(!isnan(self->cf_bond) && (self->cf_bond >= 0));
	assert(!isnan(self->cf_angle) && (self->cf_angle >= 0));
	assert(self->total_steps > 0);
	assert(self->observe_interval_mac > 0);
	assert(self->observe_interval_mic > 0);
	assert(self->root_dir != NULL);
}

int32_t getNumPtcl(const Parameter* self)
{
	return self->num_ptcl;
}

double getBondLen(const Parameter* self)
{
	return self->bond_len;
}

double getStepLen(const Parameter* self)
{
	return self->step_len;
}

int32_t getTotalSteps(const Parameter* self)
{
	return self->total_steps;
}

int32_t getObserveIntervalMac(const Parameter* self)
{
	return self->observe_interval_mac;
}

int32_t getObserveIntervalMic(const Parameter* self)
{
	return self->observe_interval_mic;
}

unsigned long getRandSeed(const Parameter* self)
{
	return self->rand_seed;
}

double getCfBond(const Parameter* self)
{
	return self->cf_bond;
}

double getCfAngle(const Parameter* self)
{
	return self->cf_angle;
}

string* getRootDir(const Parameter* self)
{
	return self->root_dir;
}

void readParameterFromFile(Parameter* self)
{
	string* input_fname = new_string_from_string(self->root_dir);
	append_char(input_fname, "/input.dat");

	// TODO: add more sophisticated way
	FILE* file = xfopen(string_to_char(input_fname), "r");
	xfscanf(file, "%d\n", &self->num_ptcl);
	xfscanf(file, "%lf\n", &self->bond_len);
	xfscanf(file, "%lf\n", &self->step_len);
	xfscanf(file, "%lf\n", &self->cf_bond);
	xfscanf(file, "%lf\n", &self->cf_angle);
	xfscanf(file, "%d\n", &self->total_steps);
	xfscanf(file, "%d\n", &self->observe_interval_mic);
	xfscanf(file, "%d\n", &self->observe_interval_mac);
	xfscanf(file, "%lu\n", &self->rand_seed);
	xfclose(file);
	delete_string(input_fname);

	checkParameterIsLoaded(self);
}