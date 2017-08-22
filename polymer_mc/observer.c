#include "observer.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "utils.h"
#include "file_utils.h"
#include "interactions.h"
#include "evolver.h"
#include "topol.h"

typedef enum {
	ENERGY = 0,
	PRESSURE,
	ACCEPT_RATIO,
	RG,
	END_TO_END,
	TRAJECT,
	BOND_LEN_DIST,

	NUM_OF_TYPES,
} ObserverType;

struct Observer_t {
	FILE* fps[NUM_OF_TYPES];
};

static const char* getFileNameFromObserverType(ObserverType type);
static void observeEnergy(Observer* self, const int32_t mc_steps, const System* system, const Boundary* bound, const Parameter* param);
static void observePressure(Observer* self, const int32_t mc_steps, const System* system, const Boundary* bound, const Parameter* param);
static void observeRg(Observer* self, const int32_t mcsteps, const System* system, const Boundary* bound, const Parameter* param);
static void observeEnd2End(Observer* self, const int32_t mcsteps, const System* system, const Boundary* bound, const Parameter* param);
static void observeAcceptRatio(Observer* self, const int32_t mcsteps, const System* system, const Parameter* param);
static void writeXYZHeader(FILE* fp, const int32_t num_ptcl, const int32_t mcsteps);
static void observeTraject(Observer* self, const int32_t mcsteps, const System* system, const Parameter* param);
static void observeBondLenDist(Observer* self, const System* system, const Boundary* bound, const Parameter* param);

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
	case TRAJECT:
		return "traject.xyz";
	case BOND_LEN_DIST:
		return "bondlen_dist.dat";
	default:
		fprintf(stderr, "ObserverType %d is not recorded.\n", type);
		exit(1);
	}
}

// TODO: implement os path join
Observer* newObserver(const string* dir_path)
{
	Observer* self = (Observer*)xmalloc(sizeof(Observer));
	for (int32_t i = 0; i < NUM_OF_TYPES; i++) {
		string* fname = new_string_from_string(dir_path);
		append_char(fname, "/");
		append_char(fname, getFileNameFromObserverType(i));
		self->fps[i] = xfopen(string_to_char(fname), "w");
	}
	return self;
}

void deleteObserver(Observer* self)
{
	for (int32_t i = 0; i < NUM_OF_TYPES; i++) {
		fclose(self->fps[i]);
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
}

#define GET_TOPOLOGY(system, param) \
	const topol* top = getTopol(system); \
	const int32_t num_bonds = getNumBonds(top); \
	const int32_t num_angles = getNumAngles(top); \
	const pair* bond_top = getBondTopol(top); \
	const triple* angle_top = getAngleTopol(top); \
	const double cf_b = getCfBond(param); \
	const double cf_a = getCfAngle(param)

static void observeEnergy(Observer* self,
	const int32_t mc_steps,
	const System* system,
	const Boundary* bound,
	const Parameter* param)
{
	static bool is_first_call = true;
	if (is_first_call) {
		fprintf(self->fps[ENERGY], "# mcsteps bond angle total \n");
		is_first_call = false;
	}

	GET_TOPOLOGY(system, param);
	const dvec* pos = getPos(system);

	double etot_bond = 0.0;
	for (int32_t b = 0; b < num_bonds; b++) {
		etot_bond += calcBondEnergy(&pos[bond_top[b].i0], &pos[bond_top[b].i1], cf_b, bound);
	}

	double etot_angle = 0.0;
	for (int32_t a = 0; a < num_angles; a++) {
		etot_angle += calcAngleEnergy(&pos[angle_top[a].i0],
			&pos[angle_top[a].i1],
			&pos[angle_top[a].i2],
			cf_a,
			bound);
	}

	const double etot = etot_bond + etot_angle;
	fprintf(self->fps[ENERGY], "%d %f %f %f\n", mc_steps, etot_bond, etot_angle, etot);
}

static void observePressure(Observer* self,
	const int32_t mcsteps,
	const System* system,
	const Boundary* bound,
	const Parameter* param)
{
	static bool is_first_call = true;
	if (is_first_call) {
		fprintf(self->fps[PRESSURE], "# mcsteps pxx   pxy   pxz   pyx   pyy   pyz   pzx   pzy   pzz\n");
		is_first_call = false;
	}

	GET_TOPOLOGY(system, param);
	const dvec* pos = getPos(system);

	dtensor3 vir_tot = {0.0};
	for (int32_t b = 0; b < num_bonds; b++) {
		const dtensor3 dvir = calcBondVirial(&pos[bond_top[b].i0], &pos[bond_top[b].i1], cf_b, bound);
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

	fprintf(self->fps[PRESSURE], "%d %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
		mcsteps,
		vir_tot.xx, vir_tot.xy, vir_tot.xz,
		vir_tot.yx, vir_tot.yy, vir_tot.yz,
		vir_tot.zx, vir_tot.zy, vir_tot.zz);
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
}

static void observeAcceptRatio(Observer* self,
	const int32_t mcsteps,
	const System* system,
	const Parameter* param)
{
	(void)param;
	fprintf(self->fps[ACCEPT_RATIO], "%d %f\n", mcsteps, getAcceptRatio(system));
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
}