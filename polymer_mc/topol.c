#include "topol.h"

#include "string_c.h"
#include "utils.h"

struct topol_t {
	int32_t num_bonds;
	pair* bond_top;
	int32_t num_angles;
	triple* angle_top;
};

static void copyPair(pair* dst, const pair* src);
static void copyTriple(triple* dst, const triple* src);
static void registerBondTopol(ptclid2topol* id2tops, const pair* bond_top);
static void registerAngleTopol(ptclid2topol* id2tops, const triple* angle_top);

topol* newTopolChain(const Parameter* param)
{
	const int32_t n = getNumPtcl(param);
	topol* top = (topol*)xmalloc(sizeof(topol));

	// bond topology
	const int32_t num_bonds = n - 1;
	top->bond_top = (pair*)xmalloc(num_bonds * sizeof(pair));
	for (int32_t i = 0; i < num_bonds; i++) {
		top->bond_top[i].i0 = i + 0;
		top->bond_top[i].i1 = i + 1;
	}
	top->num_bonds = num_bonds;

	// angle topology
	const int32_t num_angles = n - 2;
	top->angle_top = (triple*)xmalloc(num_angles * sizeof(triple));
	for (int32_t i = 0; i < num_angles; i++) {
		top->angle_top[i].i0 = i + 0;
		top->angle_top[i].i1 = i + 1;
		top->angle_top[i].i2 = i + 2;
	}
	top->num_angles = num_angles;

	return top;
}

/*topol* newTopolMesh()
{

}*/

void deleteTopol(topol* top)
{
	xfree(top->bond_top);
	xfree(top->angle_top);
	xfree(top);
}

int32_t getNumBonds(const topol* top)
{
	return top->num_bonds;
}

int32_t getNumAngles(const topol* top)
{
	return top->num_angles;
}

const pair* getBondTopol(const topol* top)
{
	return top->bond_top;
}

const triple* getAngleTopol(const topol* top)
{
	return top->angle_top;
}

static void copyPair(pair* dst, const pair* src)
{
	dst->i0 = src->i0;
	dst->i1 = src->i1;
}

static void copyTriple(triple* dst, const triple* src)
{
	dst->i0 = src->i0;
	dst->i1 = src->i1;
	dst->i2 = src->i2;
}

static void registerBondTopol(ptclid2topol* id2tops,
	const pair* bond_top)
{
	const int32_t i = bond_top->i0;
	const int32_t j = bond_top->i1;

	int32_t i_id = id2tops[i].num_pair;
	int32_t j_id = id2tops[j].num_pair;

	copyPair(&id2tops[i].pair[i_id++], bond_top);
	copyPair(&id2tops[j].pair[j_id++], bond_top);

	id2tops[i].num_pair = i_id;
	id2tops[j].num_pair = j_id;
}

static void registerAngleTopol(ptclid2topol* id2tops,
	const triple* angle_top)
{
	const int32_t i = angle_top->i0;
	const int32_t j = angle_top->i1;
	const int32_t k = angle_top->i2;

	int32_t i_id = id2tops[i].num_triple;
	int32_t j_id = id2tops[j].num_triple;
	int32_t k_id = id2tops[k].num_triple;

	copyTriple(&id2tops[i].triple[i_id++], angle_top);
	copyTriple(&id2tops[j].triple[j_id++], angle_top);
	copyTriple(&id2tops[k].triple[k_id++], angle_top);

	id2tops[i].num_triple = i_id;
	id2tops[j].num_triple = j_id;
	id2tops[k].num_triple = k_id;
}

ptclid2topol* newId2Topol(const topol* top,
	const Parameter* param)
{
	const int32_t num_ptcls = getNumPtcl(param);
	ptclid2topol* id2tops = (ptclid2topol*)xmalloc(num_ptcls * sizeof(ptclid2topol));
	for (int32_t i = 0; i < num_ptcls; i++) {
		id2tops[i].num_pair = 0;
		id2tops[i].num_triple = 0;
	}

	const int32_t num_bonds = top->num_bonds;
	for (int32_t bond = 0; bond < num_bonds; bond++) {
		registerBondTopol(id2tops, &top->bond_top[bond]);
	}

	const int32_t num_angles = top->num_angles;
	for (int32_t angle = 0; angle < num_angles; angle++) {
		registerAngleTopol(id2tops, &top->angle_top[angle]);
	}
	return id2tops;
}

void deleteId2Topol(ptclid2topol* id2tops)
{
	xfree(id2tops);
}

void debugDumpTopolInfo(const topol* top,
	const Parameter* param)
{
	const string* cdir = getRootDir(param);
	string* fname_top = new_string_from_string(cdir);
	append_char(fname_top, "/topology.dat");

	FILE* fp_top = xfopen(string_to_char(fname_top), "w");
	// bonds
	fprintf(fp_top, "# Bond topology.\n");
	const int32_t num_bonds = getNumBonds(top);
	const pair* bond_top = getBondTopol(top);
	for (int32_t i = 0; i < num_bonds; i++) {
		fprintf(fp_top, "%d %d\n",
			bond_top[i].i0, bond_top[i].i1);
	}

	// angles
	fprintf(fp_top, "# Angle topology.\n");
	const int32_t num_angles = getNumAngles(top);
	const triple* angle_top = getAngleTopol(top);
	for (int32_t i = 0; i < num_angles; i++) {
		fprintf(fp_top, "%d %d %d\n",
			angle_top[i].i0, angle_top[i].i1, angle_top[i].i2);
	}
	xfclose(fp_top);
	delete_string(fname_top);
}

void debugDumpId2TopolInfo(const ptclid2topol* id2top,
	const Parameter* param)
{
	const string* cdir = getRootDir(param);
	string* fname_id2top = new_string_from_string(cdir);
	append_char(fname_id2top, "/id2top.dat");
	FILE* fp_id2top = xfopen(string_to_char(fname_id2top), "w");
	const int32_t num_ptcls = getNumPtcl(param);
	for (int32_t i = 0; i < num_ptcls; i++) {
		fprintf(fp_id2top, "%d ", i);
		for (int32_t b = 0; b < id2top[i].num_pair; b++) {
			fprintf(fp_id2top, "(%d, %d) ",
				id2top[i].pair[b].i0, id2top[i].pair[b].i1);
		}
		for (int32_t a = 0; a < id2top[i].num_triple; a++) {
			fprintf(fp_id2top, "(%d, %d, %d) ",
				id2top[i].triple[a].i0,	id2top[i].triple[a].i1, id2top[i].triple[a].i2);
		}
		fprintf(fp_id2top, "\n");
	}
	delete_string(fname_id2top);
	xfclose(fp_id2top);
}