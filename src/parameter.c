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
  double init_blen;
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
  delete_string(self->boundary_name);
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

#define DUMP_WITH_TAG(fmt, value)               \
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
  delete_string(fname);
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

double getInitBondLen(const Parameter* self)
{
  return self->init_blen;
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

#define MATCH(value, Dtype)                                             \
  do {                                                                  \
    string* value_name = new_string_from_char(STR(value));              \
    if (eq_string(value_name, vector_ptr_string_at(keys, iter))) {      \
      STRUCT_AT(self, value) = CONCAT(string_to_, Dtype)(vector_ptr_string_at(values, iter)); \
    }                                                                   \
    delete_string(value_name);                                          \
  } while (0)

void readParameterFromFile(Parameter* self)
{
  string* input_fname = new_string_from_string(self->root_dir);
  append_char(input_fname, "/input.dat");
  FILE* fp = xfopen(string_to_char(input_fname), "r");
  delete_string(input_fname);

  vector_ptr_string* input_lines = read_lines(fp);
  vector_ptr_string* keys = vector_ptr_string_new();
  vector_ptr_string* values = vector_ptr_string_new();
  getKeysAndValues(input_lines, keys, values);

  size_t iter = 0, num_inputs = vector_ptr_string_size(keys);
  while (num_inputs != iter) {
    if (string_to_char(vector_ptr_string_at(keys, iter))[0] == ';') {
      fprintf(stdout, "# Skipping comment line.\n");
      iter++;
      continue;
    }

#ifdef SIMULATION_3D
    MATCH(side_dim_x, int32_t);
    MATCH(side_dim_y, int32_t);
    self->num_ptcl = self->side_dim_x * self->side_dim_y;
#else
    MATCH(num_ptcl, int32_t);
#endif
    MATCH(bond_len, double);
    MATCH(init_blen, double);
    MATCH(step_len, double);
    MATCH(cf_bond, double);
    MATCH(cf_angle, double);
    MATCH(total_steps, int32_t);
    MATCH(observe_interval_mic, int32_t);
    MATCH(observe_interval_mac, int32_t);
    MATCH(boundary_name, string);
    if (self->boundary_name) { /// boundary name is already set.
      if (getBoundaryTypeFromName(self->boundary_name) == PERIODIC) {
        MATCH(box_length.x, double);
        MATCH(box_length.y, double);
#ifdef SIMULATION_3D
        MATCH(box_length.z, double);
#endif
      }
    }
    MATCH(rand_seed, uint32_t);

    iter++;
  }
  delete_splitted_strings(input_lines);
  delete_splitted_strings(keys);
  delete_splitted_strings(values);
  xfclose(fp);
}
