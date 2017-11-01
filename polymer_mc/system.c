#include "system.h"

#include "topol.h"
#include "evolver.h"
#include "utils.h"
#include "file_utils.h"
#include "observer.h"
#include "mt_rand.h"

struct System_t {
  dvec* pos;
  topol* top;
  ptclid2topol* id2top;
  double accept_ratio;
};

System* newSystem(void)
{
  return (System*) xmalloc(sizeof(System));
}

void deleteSystem(System* self)
{
  xfree(self->pos);
  deleteTopol(self->top);
  deleteId2Topol(self->id2top);
  xfree(self);
}

topol* getTopol(const System* self)
{
  return self->top;
}

ptclid2topol* getPtclId2Topol(const System* self)
{
  return self->id2top;
}

dvec* getPos(const System* self)
{
  return self->pos;
}

double getAcceptRatio(const System* self)
{
  return self->accept_ratio;
}

void initializeSystem(System* self,
                      const Boundary* bound,
                      const Parameter* param,
                      confMaker conf_make,
                      topolMaker topol_make)
{
  // create initial configuration
  const int32_t num_ptcl = getNumPtcl(param);
  self->pos = (dvec*)xmalloc(num_ptcl * sizeof(dvec));
  conf_make(self, param);

  // create topology
  self->top = topol_make(param, bound);
  self->id2top = newId2Topol(self->top, param);

  debugDumpTopolInfo(self->top, param);
  debugDumpId2TopolInfo(self->id2top, param);

  // clear acceptance ratio
  self->accept_ratio = 0.0;
}

void readRestartConfig(System* self,
                       const Parameter* param)
{
  dvec* pos = getPos(self);
  const int32_t num_ptcls = getNumPtcl(param);
  const string* root_dir = getRootDir(param);
  string* fname = new_string_from_string(root_dir);
  append_char(fname, "/init_config.bin");

  int32_t n = 0;
  FILE* fp = xfopen(string_to_char(fname), "r");
  fread(&n, sizeof(int32_t), 1, fp);
  if (n != num_ptcls) {
    fprintf(stderr, "%d particles are read from init_confib.bin. However, %d is specified in input.dat\n",
            n, num_ptcls);
  }
  fread(pos, sizeof(dvec), n, fp);

  xfclose(fp);
  delete_string(fname);
}

void writeFinalConfig(System* self,
                      const Parameter* param)
{
  const int32_t num_ptcls = getNumPtcl(param);
  const dvec* pos = getPos(self);

  const string* root_dir = getRootDir(param);
  string* fname = new_string_from_string(root_dir);
  append_char(fname, "/fin_config.bin");

  FILE* fp = xfopen(string_to_char(fname), "w");
  fwrite((void *)&num_ptcls, sizeof(int32_t), 1, fp);
  fwrite((void *)pos, sizeof(dvec), num_ptcls, fp);
  xfclose(fp);

  delete_string(fname);
}

// NOTE: main simulation loop is described here.
void executeSimulation(System* self,
                       const Boundary* boundary,
                       const Parameter* param)
{
  Observer* observer = newObserver(getRootDir(param));
  MTstate* mtst = newMTstate();
  init_genrand(mtst, getRandSeed(param));

  const int32_t tot_steps = getTotalSteps(param);
  const int32_t observe_interval_mic = getObserveIntervalMic(param);
  const int32_t observe_interval_mac = getObserveIntervalMac(param);
  for (int32_t i = 0; i < tot_steps; i++) {
    self->accept_ratio = evolveMc(self, param, boundary, mtst);
    if (i % observe_interval_mic == 0) observeMicroVars(observer, i, self, boundary, param);
    if (i % observe_interval_mac == 0) observeMacroVars(observer, i, self, boundary, param);
  }

  deleteObserver(observer);
  deleteMTstate(mtst);
}
