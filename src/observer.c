#include "observer.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>

#include "utils.h"
#include "file_utils.h"
#include "interactions.h"
#include "evolver.h"
#include "topol.h"
#include "math_utils.h"

typedef enum {
  ENERGY = 0,
  PRESSURE,
  ACCEPT_RATIO,
  RG,
  END_TO_END,
  BOND_LEN_DIST,
  TRAJECT,
  FLUCT_SPECTRUM,

  NUM_OF_TYPES,
} ObserverType;

typedef void (* finalizer_t)(Observer* self);

struct Observer_t {
  FILE* fps[NUM_OF_TYPES];
  int32_t num_frames[NUM_OF_TYPES];
  void* buffer[NUM_OF_TYPES];
  finalizer_t finalizer[NUM_OF_TYPES];
};

static const char* getFileNameFromObserverType(ObserverType type);

static void initializeEnergyObserver(Observer* self);
static void finalizeEnergyObserver(Observer* self);
static void observeEnergy(Observer* self, const int32_t mc_steps, const System* system, const Boundary* bound, const Parameter* param);

static void initializePressureObserver(Observer* self);
static void finalizePressureObserver(Observer* self);
static void observePressure(Observer* self, const int32_t mc_steps, const System* system, const Boundary* bound, const Parameter* param);

static void observeRg(Observer* self, const int32_t mcsteps, const System* system, const Boundary* bound, const Parameter* param);
static void observeEnd2End(Observer* self, const int32_t mcsteps, const System* system, const Boundary* bound, const Parameter* param);
static void observeAcceptRatio(Observer* self, const int32_t mcsteps, const System* system, const Parameter* param);
static void writeXYZHeader(FILE* fp, const int32_t num_ptcl, const int32_t mcsteps);
static void observeTraject(Observer* self, const int32_t mcsteps, const System* system, const Parameter* param);
static void observeBondLenDist(Observer* self, const System* system, const Boundary* bound, const Parameter* param);

struct SpectrumBuffer_t;
typedef struct SpectrumBuffer_t SpectrumBuffer;
static double complex calcftilde1D(const dvec* pos, const int32_t n_ptcl, const double qx_at);
static double complex calcftilde2D(const dvec* pos, const int32_t n_ptcl,
                                   const double qx_at, const double qy_at);
static void doFourierTransform1D(const dvec* pos, const int32_t n_ptcl, SpectrumBuffer* sbuffer);
static void doFourierTransform2D(const dvec* pos, const int32_t n_ptcl, SpectrumBuffer* sbuffer);
static void setQvector1D(double* qx, const double q_low, const double q_up, const int32_t ndiv);
static void setQvector2D(double* qx, double* qy,
                       const double qx_low, const double qx_up,
                       const double qy_low, const double qy_up,
                       const int32_t nx_div, const int32_t ny_div);
static void initializeFluctSpetrumObserver(Observer* self, const Parameter* param);
static void finalizeFluctSpetrumObserver(Observer* self);
static void observeFluctSpectrum(Observer* self, const System* system, const Parameter* param);

static const char* getFileNameFromObserverType(ObserverType type)
{
  switch (type) {
  case ENERGY:
    return "energy.dat";
  case PRESSURE:
    return "pressure.dat";
  case ACCEPT_RATIO:
    return "accept_ratio.dat";
  case RG:
    return "rg.dat";
  case END_TO_END:
    return "end2end.dat";
  case BOND_LEN_DIST:
    return "bondlen_dist.dat";
  case TRAJECT:
    return "traject.xyz";
  case FLUCT_SPECTRUM:
    return "fluct_spectrum.dat";
  default:
    fprintf(stderr, "ObserverType %d is not recorded.\n", type);
    fprintf(stderr, "Error occurs at %s %d\n", __FILE__, __LINE__);
    exit(1);
  }
}

// TODO: implement os path join
Observer* newObserver(const string* dir_path)
{
  Observer* self = (Observer*)xmalloc(sizeof(Observer));
  string* fname = new_string();
  for (int32_t i = 0; i < NUM_OF_TYPES; i++) {
    append_string(fname, dir_path);
    append_char(fname, "/");
    append_char(fname, getFileNameFromObserverType(i));
    self->fps[i]         = xfopen(string_to_char(fname), "w");
    self->num_frames[i]  = 0;
    self->buffer[i]      = NULL;
    self->finalizer[i]   = NULL;
    clear_string(fname);
  }
  delete_string(fname);
  return self;
}

void deleteObserver(Observer* self)
{
  for (int32_t i = 0; i < NUM_OF_TYPES; i++) {
    if (self->finalizer[i]) {
      self->finalizer[i](self);
    }
  }
  for (int32_t i = 0; i < NUM_OF_TYPES; i++) {
    fclose(self->fps[i]);
  }
  for (int32_t i = 0; i < NUM_OF_TYPES; i++) {
    xfree(self->buffer[i]);
  }
  xfree(self);
}

void observeMicroVars(Observer* self,
                      const int32_t mc_steps,
                      const System* system,
                      const Boundary* bound,
                      const Parameter* param)
{
  observeTraject(self, mc_steps, system, param);
  observeBondLenDist(self, system, bound, param);
}

void observeMacroVars(Observer* self,
                      const int32_t mc_steps,
                      const System* system,
                      const Boundary* bound,
                      const Parameter* param)
{
  observeEnergy(self, mc_steps, system, bound, param);
  observePressure(self, mc_steps, system, bound, param);
  observeRg(self, mc_steps, system, bound, param);
  observeEnd2End(self, mc_steps, system, bound, param);
  observeAcceptRatio(self, mc_steps, system, param);
  if (getBoundaryType(bound) == PERIODIC) {
    observeFluctSpectrum(self, system, param);
  }
}

#define GET_TOPOLOGY(system, param)             \
  const topol* top = getTopol(system);          \
  const int32_t num_bonds = getNumBonds(top);   \
  const int32_t num_angles = getNumAngles(top); \
  const pair* bond_top = getBondTopol(top);     \
  const triple* angle_top = getAngleTopol(top); \
  const double cf_b = getCfBond(param);         \
  const double cf_a = getCfAngle(param);        \
  const double l0  = getBondLen(param)

typedef struct EnergyBuffer_t {
  double bond;
  double angle;
  double total;
} EnergyBuffer;

static void initializeEnergyObserver(Observer* self)
{
  fprintf(self->fps[ENERGY], "# mcsteps bond angle total \n");
  self->buffer[ENERGY] = (EnergyBuffer*) xmalloc(sizeof(EnergyBuffer));
  EnergyBuffer* ebuffer = (EnergyBuffer*) self->buffer[ENERGY];
  ebuffer->bond = ebuffer->angle = ebuffer->total = 0.0;
  self->finalizer[ENERGY] = finalizeEnergyObserver;
}

static void finalizeEnergyObserver(Observer* self)
{
  EnergyBuffer* ebuffer = (EnergyBuffer*) self->buffer[ENERGY];
  ebuffer->bond  /= self->num_frames[ENERGY];
  ebuffer->angle /= self->num_frames[ENERGY];
  ebuffer->total /= self->num_frames[ENERGY];
  fprintf(self->fps[ENERGY], "# mean = %f %f %f\n",
          ebuffer->bond, ebuffer->angle, ebuffer->total);
}

static void observeEnergy(Observer* self,
                          const int32_t mc_steps,
                          const System* system,
                          const Boundary* bound,
                          const Parameter* param)
{
  static bool is_first_call = true;
  if (is_first_call) {
    initializeEnergyObserver(self);
    is_first_call = false;
  }

  GET_TOPOLOGY(system, param);
  const dvec* pos = getPos(system);

  // sum bonded energy
  double etot_bond = 0.0;
  for (int32_t b = 0; b < num_bonds; b++) {
    etot_bond += calcBondEnergy(&pos[bond_top[b].i0], &pos[bond_top[b].i1], cf_b, l0, bound);
  }

  // sum angle energy
  double etot_angle = 0.0;
  for (int32_t a = 0; a < num_angles; a++) {
    etot_angle += calcAngleEnergy(&pos[angle_top[a].i0],
                                  &pos[angle_top[a].i1],
                                  &pos[angle_top[a].i2],
                                  cf_a,
                                  bound);
  }

  // print out energy
  const double etot = etot_bond + etot_angle;
  fprintf(self->fps[ENERGY], "%d %f %f %f\n", mc_steps, etot_bond, etot_angle, etot);

  // accumulate result
  EnergyBuffer* ebuffer = (EnergyBuffer*) self->buffer[ENERGY];
  ebuffer->bond += etot_bond;
  ebuffer->angle += etot_angle;
  ebuffer->total += etot;
  self->num_frames[ENERGY]++;
}

typedef struct PressureBuffer_t {
  dtensor3 virial;
} PressureBuffer;

static void initializePressureObserver(Observer* self)
{
  self->buffer[PRESSURE] = (PressureBuffer*) xmalloc(sizeof(PressureBuffer));
  PressureBuffer* pbuffer = (PressureBuffer*) self->buffer[PRESSURE];
  dtensor3_clear(&pbuffer->virial);
  self->finalizer[PRESSURE] = finalizePressureObserver;
}

static void finalizePressureObserver(Observer* self)
{
  PressureBuffer* pbuffer = (PressureBuffer*) self->buffer[PRESSURE];
  dtensor3_mul_scalar(&pbuffer->virial, 1.0 / self->num_frames[PRESSURE]);
  fprintf(self->fps[PRESSURE],
          "# mean = %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
          pbuffer->virial.xx, pbuffer->virial.xy, pbuffer->virial.xz,
          pbuffer->virial.yx, pbuffer->virial.yy, pbuffer->virial.yz,
          pbuffer->virial.zx, pbuffer->virial.zy, pbuffer->virial.zz);
}

static void observePressure(Observer* self,
                            const int32_t mcsteps,
                            const System* system,
                            const Boundary* bound,
                            const Parameter* param)
{
  static bool is_first_call = true;
  if (is_first_call) {
    fprintf(self->fps[PRESSURE],
            "# mcsteps pxx   pxy   pxz   pyx   pyy   pyz   pzx   pzy   pzz\n");
    initializePressureObserver(self);
    is_first_call = false;
  }

  GET_TOPOLOGY(system, param);
  const dvec* pos = getPos(system);

  dtensor3 vir_tot = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  for (int32_t b = 0; b < num_bonds; b++) {
    const dtensor3 dvir = calcBondVirial(&pos[bond_top[b].i0],
                                         &pos[bond_top[b].i1],
                                         cf_b,
                                         l0,
                                         bound);
    dtensor3_add(&vir_tot, &dvir);
  }
  for (int32_t a = 0; a < num_angles; a++) {
    const dtensor3 dvir = calcAngleVirial(&pos[angle_top[a].i0],
                                          &pos[angle_top[a].i1],
                                          &pos[angle_top[a].i2],
                                          cf_a,
                                          bound);
    dtensor3_add(&vir_tot, &dvir);
  }
  fprintf(self->fps[PRESSURE],
          "%d %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
          mcsteps,
          vir_tot.xx, vir_tot.xy, vir_tot.xz,
          vir_tot.yx, vir_tot.yy, vir_tot.yz,
          vir_tot.zx, vir_tot.zy, vir_tot.zz);

  // accumulate results
  PressureBuffer* pbuffer = (PressureBuffer*) self->buffer[PRESSURE];
  dtensor3_add(&pbuffer->virial, &vir_tot);
  self->num_frames[PRESSURE]++;
}

static void observeRg(Observer* self,
                      const int32_t mcsteps,
                      const System* system,
                      const Boundary* bound,
                      const Parameter* param)
{
  const dvec* pos = getPos(system);
  const int32_t num_ptcl = getNumPtcl(param);

  dvec cmpos = {.x = 0.0, .y = 0.0, .z = 0.0};
  for (int32_t i = 0; i < num_ptcl; i++) {
    add_dvec(&cmpos, &pos[i]);
  }
  div_scalar(&cmpos, (double)num_ptcl);

  double rg = 0.0;
  for (int32_t i = 0; i < num_ptcl; i++) {
    rg += distance2(&cmpos, &pos[i], bound);
  }
  rg /= num_ptcl;
  rg = sqrt(rg);

  fprintf(self->fps[RG], "%d %f\n", mcsteps, rg);
  self->num_frames[RG]++;
}

static void observeEnd2End(Observer* self,
                           const int32_t mcsteps,
                           const System* system,
                           const Boundary* bound,
                           const Parameter* param)
{
  const dvec* pos = getPos(system);
  const int32_t num_ptcl = getNumPtcl(param);
  const double e2e = distance(&pos[0], &pos[num_ptcl - 1], bound);
  fprintf(self->fps[END_TO_END], "%d %f\n", mcsteps, e2e);
  self->num_frames[END_TO_END]++;
}

static void observeAcceptRatio(Observer* self,
                               const int32_t mcsteps,
                               const System* system,
                               const Parameter* param)
{
  UNUSED_PARAMETER(param);
  fprintf(self->fps[ACCEPT_RATIO], "%d %f\n", mcsteps, getAcceptRatio(system));
  self->num_frames[ACCEPT_RATIO]++;
}

static void writeXYZHeader(FILE* fp,
                           const int32_t num_ptcl,
                           const int32_t mcsteps)
{
  fprintf(fp, "%d\n", num_ptcl);
  fprintf(fp, "mcsteps = %d\n", mcsteps);
}

static void observeTraject(Observer* self,
                           const int32_t mcsteps,
                           const System* system,
                           const Parameter* param)
{
  const dvec* pos = getPos(system);
  const int32_t num_ptcl = getNumPtcl(param);
  writeXYZHeader(self->fps[TRAJECT], num_ptcl, mcsteps);
  for (int32_t i = 0; i < num_ptcl; i++) {
    fprintf(self->fps[TRAJECT],
            "C %.15g %.15g %.15g\n", pos[i].x, pos[i].y, pos[i].z);
  }
  self->num_frames[TRAJECT]++;
}

static void observeBondLenDist(Observer* self,
                               const System* system,
                               const Boundary* bound,
                               const Parameter* param)
{
  const dvec* pos = getPos(system);
  const int32_t num_ptcl = getNumPtcl(param);
  for (int32_t i = 1; i < num_ptcl; i++) {
    fprintf(self->fps[BOND_LEN_DIST],
            "%.15g\n", distance(&pos[i], &pos[i-1], bound));
  }
  self->num_frames[BOND_LEN_DIST]++;
}

typedef struct SpectrumBuffer_t {
  double *qx, *qy;
  int32_t nx_div, ny_div;
  double complex* spect_sum;
  double factor; // normalize factor
} SpectrumBuffer;

static double complex calcftilde1D(const dvec* pos,
                                   const int32_t n_ptcl,
                                   const double qx_at)
{
  double complex sum = 0.0 + 0.0 * I;
  for (int32_t i = 0; i < n_ptcl; i++) {
    sum += pos[i].y * cexp(-I * qx_at * pos[i].x);
  }
  return sum;
}

static double complex calcftilde2D(const dvec* pos,
                                   const int32_t n_ptcl,
                                   const double qx_at,
                                   const double qy_at)
{
  double complex sum = 0.0 + 0.0 * I;
  for (int32_t i = 0; i < n_ptcl; i++) {
    sum += pos[i].z * cexp(-I * (qx_at * pos[i].x + qy_at * pos[i].y));
  }
  return sum;
}

static void doFourierTransform1D(const dvec* pos,
                                 const int32_t n_ptcl,
                                 SpectrumBuffer* sbuffer)
{
  const int32_t nx_div = sbuffer->nx_div;
  for (int32_t i = 0; i < nx_div; i++) {
    sbuffer->spect_sum[i] += calcftilde1D(pos, n_ptcl, sbuffer->qx[i]);
  }
}

static void doFourierTransform2D(const dvec* pos,
                                 const int32_t n_ptcl,
                                 SpectrumBuffer* sbuffer)
{
  const int32_t nx_div = sbuffer->nx_div;
  const int32_t ny_div = sbuffer->ny_div;
  int32_t cnt = 0;
  for (int32_t iy = 0; iy < ny_div; iy++) {
    for (int32_t ix = 0; ix < nx_div; ix++) {
      sbuffer->spect_sum[cnt++] += calcftilde2D(pos, n_ptcl,
                                                sbuffer->qx[ix],
                                                sbuffer->qy[iy]);
    }
  }
}

static void setQvector1D(double* qx,
                       const double q_low,
                       const double q_up,
                       const int32_t ndiv)
{
  const double dq = (q_up - q_low) / ndiv;
  for (int32_t i = 0; i < ndiv; i++) {
    qx[i] = (i + 0.5) * dq;
  }
}

static void setQvector2D(double* qx, double* qy,
                       const double qx_low, const double qx_up,
                       const double qy_low, const double qy_up,
                       const int32_t nx_div, const int32_t ny_div)
{
  const double dqx = (qx_up - qx_low) / nx_div;
  const double dqy = (qy_up - qy_low) / ny_div;
  for (int32_t iy = 0; iy < ny_div; iy++) {
    for (int32_t ix = 0; ix < nx_div; ix++) {
      qx[ix] = (ix + 0.5) * dqx;
      qy[iy] = (iy + 0.5) * dqy;
    }
  }
}

static void initializeFluctSpetrumObserver(Observer* self,
                                           const Parameter* param)
{
  self->buffer[FLUCT_SPECTRUM] = (SpectrumBuffer*) xmalloc(sizeof(SpectrumBuffer));
  SpectrumBuffer* sbuffer      = (SpectrumBuffer*) self->buffer[FLUCT_SPECTRUM];
  const dvec box_length        = getBoxlength(param);
  const double b_len = getBondLen(param);
#ifdef SIMULATION_3D
  sbuffer->nx_div = 50;
  sbuffer->ny_div = 50;
  sbuffer->factor = 1.0 / sqrt(box_length.x * box_length.y);
  sbuffer->qx = (double*) xmalloc(sbuffer->nx_div * sizeof(double));
  sbuffer->qy = (double*) xmalloc(sbuffer->ny_div * sizeof(double));
#else
  sbuffer->nx_div = 100;
  sbuffer->ny_div = 1;
  sbuffer->factor = 1.0 / sqrt(box_length.x);
  sbuffer->qx = (double*) xmalloc(sbuffer->nx_div * sizeof(double));
  sbuffer->qy = NULL;
#endif
  const int32_t ndiv = sbuffer->nx_div * sbuffer->ny_div;
  sbuffer->spect_sum
    = (double complex*) xmalloc(ndiv * sizeof(double complex));

#ifdef SIMULATION_3D
  const double qx_low = 1.0 / (M_PI * box_length.x);
  const double qy_low = 1.0 / (M_PI * box_length.y);
  const double qx_up  = 1.0 / (M_PI * b_len), qy_up  = 1.0 / (M_PI * b_len);
  setQvector2D(sbuffer->qx, sbuffer->qy,
               qx_low, qx_up,
               qy_low, qy_up,
               sbuffer->nx_div, sbuffer->ny_div);
  UNUSED_FUNCTION(setQvector1D);
#else
  const double qx_low = 1.0 / (M_PI * box_length.x);
  const double qx_up  = 1.0 / (M_PI * b_len);
  setQvector1D(sbuffer->qx,
               qx_low, qx_up,
               sbuffer->nx_div);
  UNUSED_FUNCTION(setQvector2D);
#endif

  for (int32_t i = 0; i < ndiv; i++) {
    sbuffer->spect_sum[i] = 0.0 + 0.0 * I;
  }
  self->finalizer[FLUCT_SPECTRUM] = finalizeFluctSpetrumObserver;
}

static void finalizeFluctSpetrumObserver(Observer* self)
{
  SpectrumBuffer* sbuffer = (SpectrumBuffer*) self->buffer[FLUCT_SPECTRUM];
  const double cf = sbuffer->factor / self->num_frames[FLUCT_SPECTRUM];
#ifdef SIMULATION_3D
  int32_t cnt = 0;
  for (int32_t iy = 0; iy < sbuffer->ny_div; iy++) {
    for (int32_t ix = 0; ix < sbuffer->nx_div; ix++) {
      const double q_norm = sqrt(sbuffer->qx[i] * sbuffer->qx[i]
                                 + sbuffer->qy[i] * sbuffer->qy[i]);
      sbuffer->spect_sum[cnt] *= cf;
      const double spect_norm
        = creal(sbuffer->spect_sum[cnt] * conj(sbuffer->spect_sum[cnt]));
      fprintf(self->fps[FLUCT_SPECTRUM],
              "%.10g %.10g\n",
              q_norm, spect_norm);
      cnt++;
    }
  }
#else
  for (int32_t i = 0; i < sbuffer->nx_div; i++) {
    sbuffer->spect_sum[i] *= cf;
    const double spect_norm
      = creal(sbuffer->spect_sum[i] * conj(sbuffer->spect_sum[i]));
    fprintf(self->fps[FLUCT_SPECTRUM],
            "%.10g %.10g\n",
            sbuffer->qx[i], spect_norm);
  }
#endif
  xfree(sbuffer->qx);
  xfree(sbuffer->qy);
  xfree(sbuffer->spect_sum);
}

static void observeFluctSpectrum(Observer* self,
                                 const System* system,
                                 const Parameter* param)
{
  static bool is_first_call = true;
  if (is_first_call) {
    initializeFluctSpetrumObserver(self, param);
    is_first_call = false;
  }

  const dvec* pos = getPos(system);
  const int32_t num_ptcl = getNumPtcl(param);
  SpectrumBuffer* sbuffer = (SpectrumBuffer*) self->buffer[FLUCT_SPECTRUM];

#ifdef SIMULATION_3D
  doFourierTransform2D(pos, num_ptcl, sbuffer);
  UNUSED_FUNCTION(doFourierTransform1D);
#else
  doFourierTransform1D(pos, num_ptcl, sbuffer);
  UNUSED_FUNCTION(doFourierTransform2D);
#endif

  self->num_frames[FLUCT_SPECTRUM]++;
}
