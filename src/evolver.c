#include "evolver.h"

#if _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "mt_rand.h"
#include "parameter.h"
#include "system.h"
#include "topol.h"
#include "utils.h"
#include "interactions.h"
#include "boundary.h"
#include "vector3.h"

static dvec kickParticle(const dvec *pos0,
                         const double disp,
                         MTstate *mtst,
                         const Boundary *bound)
{
  dvec disp_vec;
  disp_vec.x = disp * (2.0 * genrand_res53(mtst) - 1.0);
  disp_vec.y = disp * (2.0 * genrand_res53(mtst) - 1.0);
#ifdef SIMULATION_3D
  disp_vec.z = disp * (2.0 * genrand_res53(mtst) - 1.0);
#else
  disp_vec.z = 0.0;
#endif
  dvec new_pos = add_dvec_new(pos0, &disp_vec);
  applyBoundaryCond(bound, &new_pos);
  return new_pos;
}

static bool newStateIsAccepted(const double deltaE,
                               MTstate *mtst)
{
  if (deltaE < 0.0)
  {
    return true;
  }
  else
  {
    const double uni_rand = genrand_res53(mtst);
    const double threshold = exp(-deltaE);
    return (uni_rand < threshold);
  }
}

static double calcBondEnergyLocalSum(dvec *pos,
                                     const int32_t id_picked,
                                     const ptclid2topol *id2top,
                                     const Boundary *bound,
                                     const double cf_bond,
                                     const double l0)
{
  double esum = 0.0;
  const int32_t num_bonds = id2top[id_picked].num_pair;
  for (int32_t bond = 0; bond < num_bonds; bond++)
  {
    const int32_t i = id2top[id_picked].pair[bond].i0;
    const int32_t j = id2top[id_picked].pair[bond].i1;
    esum += calcBondEnergy(&pos[i], &pos[j], cf_bond, l0, bound);
  }
  return esum;
}

static double calcAngleEnergyLocalSum(dvec *pos,
                                      const int32_t id_picked,
                                      const ptclid2topol *id2top,
                                      const Boundary *bound,
                                      const double cf_angle)
{
  double esum = 0.0;
  const int32_t num_angles = id2top[id_picked].num_triple;
  for (int32_t angle = 0; angle < num_angles; angle++)
  {
    const int32_t i = id2top[id_picked].triple[angle].i0;
    const int32_t j = id2top[id_picked].triple[angle].i1;
    const int32_t k = id2top[id_picked].triple[angle].i2;
    esum += calcAngleEnergy(&pos[i], &pos[j], &pos[k], cf_angle, bound);
  }
  return esum;
}

static double calcLocEnergy(dvec *pos,
                            const int32_t id_picked,
                            const ptclid2topol *id2top,
                            const Boundary *bound,
                            const double cf_bond,
                            const double cf_angle,
                            const double l0)
{
  return calcBondEnergyLocalSum(pos, id_picked, id2top, bound, cf_bond, l0) + calcAngleEnergyLocalSum(pos, id_picked, id2top, bound, cf_angle);
}

static void mcStep(dvec *pos,
                   MTstate *mtst,
                   int32_t *num_accepted,
                   const ptclid2topol *id2top,
                   const Boundary *bound,
                   const double disp,
                   const double cf_bond,
                   const double cf_angle,
                   const double l0,
                   const int32_t id_lo,
                   const int32_t id_hi)
{
  const int32_t id_picked = genrand_int31_range(mtst, id_lo, id_hi);
  const dvec pos_tmp = pos[id_picked];

  const double e_locsum_bef = calcLocEnergy(pos, id_picked, id2top, bound, cf_bond, cf_angle, l0);
  pos[id_picked] = kickParticle(&pos_tmp, disp, mtst, bound);
  const double e_locsum_aft = calcLocEnergy(pos, id_picked, id2top, bound, cf_bond, cf_angle, l0);
  const double dE = e_locsum_aft - e_locsum_bef;

  if (newStateIsAccepted(dE, mtst))
  {
    (*num_accepted)++;
  }
  else
  {
    pos[id_picked] = pos_tmp;
  }
}

double evolveMc(System *system,
                const Parameter *param,
                const Boundary *bound,
                MTstate *mtst)
{
  const int32_t num_ptcl = getNumPtcl(param);
  const double step_len = getStepLen(param);
  const double cf_bond = getCfBond(param);
  const double cf_angle = getCfAngle(param);
  const double l0 = getBondLen(param);

  for (int32_t p = 0; p < num_ptcl; p++)
  {
    const double PIVOT_MOVE_PROBABILITY = 0.1; // let's roll with 10% for now
    double rand_val = genrand_real1(mtst);
    if (rand_val < PIVOT_MOVE_PROBABILITY)
    {
      continue;
    }

    while (1)
    {
      int32_t pivot = 1 + genrand_int31(mtst) % (num_ptcl - 2);
      // hm. lol.
      double axis[3] = {genrand_real1(mtst), genrand_real1(mtst), genrand_real1(mtst)};
      double angle = 2 * M_PI * genrand_real1(mtst); // Full rotation range

      for (int32_t i = pivot + 1; i < num_ptcl; i++)
      {
        printf("rotate the particle %d around the axis (%f, %f, %f) by %f\n", i, axis[0], axis[1], axis[2], angle);
      }
    }

    /*     mcStep(getPos(system), mtst, NULL, getPtclId2Topol(system), bound,
              step_len, cf_bond, cf_angle, l0,
               0, num_ptcl - 1); */
  }

  int32_t id_movable_lo = 0, id_movable_hi = num_ptcl - 1;
  if (getBoundaryType(bound) == PERIODIC)
  {
    id_movable_lo++;
    id_movable_hi--;
  }

  dvec *pos = getPos(system);
  ptclid2topol *id2top = getPtclId2Topol(system);

  int32_t num_accepted = 0;
  for (int32_t p = 0; p < num_ptcl; p++)
  {
    mcStep(pos, mtst, &num_accepted, id2top, bound,
           step_len, cf_bond, cf_angle, l0,
           id_movable_lo, id_movable_hi);
  }

  return (double)num_accepted / (double)num_ptcl;
}

// TODO(erwan): this seems directionally correct, but prob not quite right
bool legalParticule(const dvec *pos, int32_t num_ptcl, const Boundary *bound)
{
  const double MIN_DISTANCE_MONOMER = 1.0; // really a placeholder for now

  for (int32_t i = 0; i < num_ptcl - 1; i++)
  {
    for (int32_t j = i + 1; j < num_ptcl; j++)
    {
      if (checkParticleOverlap(pos, num_ptcl, bound))
      {
        return false;
      }
    }
  }

  if (bound == PERIODIC)
  {
    for (int32_t i = 0; i < num_ptcl; i++)
    {
      if (!particuleBoundary(pos[i], bound))
      {
        return false;
      }
    }
  }

  return true;
}
bool particuleBoundary(dvec point, const Boundary *bound)
{
  if (getBoundaryType(bound) == FREE)
  {
    return true;
  }
  else if (getBoundaryType(bound) == PERIODIC)
  {
    dvec boxLength;

    return (point.x >= 0 && point.x <= boxLength.x &&
            point.y >= 0 && point.y <= boxLength.y &&
            point.z >= 0 && point.z <= boxLength.z);
  }

  return false;
}
double boundaryDistance(const dvec pos1, const dvec pos2, const Boundary *bound)
{
  return distance(&pos1, &pos2, bound); 
}

bool checkParticleOverlap(const dvec *pos, int32_t num_ptcl, const Boundary *bound)
{
  const double MIN_DISTANCE = 1.0; // this is a placeholder for now

  for (int32_t i = 0; i < num_ptcl - 1; i++)
  {
    for (int32_t j = i + 1; j < num_ptcl; j++)
    {
      if (boundaryDistance(pos[i], pos[j], bound) < MIN_DISTANCE)
      {
        return true;
      }
    }
  }
  return false;
}
