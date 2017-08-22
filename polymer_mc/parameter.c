#include "parameter.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "string_c.h"
#include "file_utils.h"
#include "utils.h"
#include "boundary.h"

struct Parameter_t {
	string* root_dir;
	int32_t num_ptcl;
	int32_t side_dim_x;
	int32_t side_dim_y;
	int32_t total_steps;
	int32_t observe_interval_mic;
	int32_t observe_interval_mac;
	double bond_len;
	double step_len;
	double cf_bond;
	double cf_angle;
	dvec box_length;
	string* boundary_name;
	uint32_t rand_seed;
};

static void initializeParameter(Parameter* self);
static void dumpAllParameter(Parameter* self);

Parameter* newParameter(const char* path)
{
	Parameter* self = (Parameter*)xmalloc(sizeof(Parameter));
	initializeParameter(self);
	self->root_dir = new_string_from_char(path);
	return self;
}

void deleteParameter(Parameter* self)
{
	dumpAllParameter(self);
	delete_string(self->root_dir);
	xfree(self);
}

static void initializeParameter(Parameter* self)
{
	self->root_dir = NULL;
	self->num_ptcl = -1;
	self->side_dim_x = -1;
	self->side_dim_y = -1;
	self->total_steps = -1;
	self->observe_interval_mic = -1;
	self->observe_interval_mac = -1;
	self->bond_len = nan("");
	self->step_len = nan("");
	self->cf_bond = nan("");
	self->cf_angle = nan("");
	self->box_length.x = nan("");
	self->box_length.y = nan("");
	self->box_length.z = 0.0;
	self->rand_seed = 0xffffffff;
	self->boundary_name = NULL;
}

#define DUMP_WITH_TAG(fmt, value) \
fprintf(fp, fmt, STR(value), self->value);

static void dumpAllParameter(Parameter* self)
{
	string* fname = new_string_from_string(self->root_dir);
	append_char(fname, "/all_param.dat");
	FILE* fp = xfopen(string_to_char(fname), "w");
	fprintf(fp, "%s = %s\n", "root_dir", string_to_char(self->root_dir));
	DUMP_WITH_TAG("%s = %d\n", num_ptcl);
	DUMP_WITH_TAG("%s = %d\n", side_dim_x);
	DUMP_WITH_TAG("%s = %d\n", side_dim_y);
	DUMP_WITH_TAG("%s = %d\n", total_steps);
	DUMP_WITH_TAG("%s = %d\n", observe_interval_mac);
	DUMP_WITH_TAG("%s = %d\n", observe_interval_mic);
	DUMP_WITH_TAG("%s = %lf\n", bond_len);
	DUMP_WITH_TAG("%s = %lf\n", step_len);
	DUMP_WITH_TAG("%s = %lf\n", cf_bond);
	DUMP_WITH_TAG("%s = %lf\n", cf_angle);
	DUMP_WITH_TAG("%s = %lf\n", box_length.x);
	DUMP_WITH_TAG("%s = %lf\n", box_length.y);
	DUMP_WITH_TAG("%s = %lf\n", box_length.z);
	DUMP_WITH_TAG("%s = %u\n", rand_seed);
	fprintf(fp, "%s = %s\n", "boundary_name", string_to_char(self->boundary_name));
	xfclose(fp);
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

uint32_t getRandSeed(const Parameter* self)
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

const string* getRootDir(const Parameter* self)
{
	return self->root_dir;
}

const string* getBoundaryName(const Parameter* self)
{
	return self->boundary_name;
}

dvec getBoxlength(const Parameter* self)
{
	return self->box_length;
}

int32_t getSideDimx(const Parameter* self)
{
	return self->side_dim_x;
}

int32_t getSideDimy(const Parameter* self)
{
	return self->side_dim_y;
}

static void getKeysAndValues(const vector_ptr_string* key_values,
	vector_ptr_string* keys,
	vector_ptr_string* values)
{
	const size_t num_key_values = vector_ptr_string_size(key_values);
	for (size_t i = 0; i < num_key_values; i++) {
		vector_ptr_string* key_value
			= split_string(vector_ptr_string_at(key_values, i), " ");
		vector_ptr_string_push_back(keys, vector_ptr_string_at(key_value, 0));
		vector_ptr_string_push_back(values, vector_ptr_string_at(key_value, 1));
		vector_ptr_string_delete(key_value);
	}
}

static size_t findElementInVectorString(const vector_ptr_string* array,
	const char* elem)
{
	const string* elem_string = new_string_from_char(elem);
	const size_t size_array = vector_ptr_string_size(array);
	for (size_t i = 0; i < size_array; i++) {
		if (eq_string(elem_string, vector_ptr_string_at(array, i))) {
			return i; // found!
		}
	}
	return size_array; // not found!
}

#define MATCH(keys, values, value_name, param, Dtype)\
do {\
	const size_t size_keys = vector_ptr_string_size(keys);\
	const size_t value_at  = findElementInVectorString(keys, STR(value_name));\
	if (value_at != size_keys) {\
		STRUCT_AT(param, value_name) = CONCAT(string_to_, Dtype)(vector_ptr_string_at(values, value_at));\
	} else {\
		fprintf(stderr, "Cannot find tag %s\n", STR(value_name));\
		exit(EXIT_FAILURE);\
	}\
} while (0)

void readParameterFromFile(Parameter* self,
	const bool is3d)
{
	string* input_fname = new_string_from_string(self->root_dir);
	append_char(input_fname, "/input.dat");
	FILE* fp = xfopen(string_to_char(input_fname), "r");
	delete_string(input_fname);

	vector_ptr_string* input_lines = read_lines(fp);
	vector_ptr_string* keys = vector_ptr_string_new();
	vector_ptr_string* values = vector_ptr_string_new();
	getKeysAndValues(input_lines, keys, values);

	if (is3d) {
		MATCH(keys, values, side_dim_x, self, int32_t);
		MATCH(keys, values, side_dim_y, self, int32_t);
		self->num_ptcl = self->side_dim_x * self->side_dim_y;
	} else {
		MATCH(keys, values, num_ptcl, self, int32_t);
	}
	MATCH(keys, values, bond_len, self, double);
	MATCH(keys, values, step_len, self, double);
	MATCH(keys, values, cf_bond, self, double);
	MATCH(keys, values, cf_angle, self, double);
	MATCH(keys, values, total_steps, self, int32_t);
	MATCH(keys, values, observe_interval_mic, self, int32_t);
	MATCH(keys, values, observe_interval_mac, self, int32_t);
	MATCH(keys, values, boundary_name, self, string);
	if (getBoundaryTypeFromName(self->boundary_name) == PERIODIC) {
		MATCH(keys, values, box_length.x, self, double);
		MATCH(keys, values, box_length.y, self, double);
		if (is3d) { MATCH(keys, values, box_length.z, self, double); }
	}
	MATCH(keys, values, rand_seed, self, uint32_t);

	delete_splitted_strings(input_lines);
	delete_splitted_strings(keys);
	delete_splitted_strings(values);

	xfclose(fp);
}