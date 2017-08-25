#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "user_defs.h"
#include "system.h"
#include "topol.h"
#include "parameter.h"
#include "config_maker.h"
#include "boundary.h"

static void check_args(const int argc,
	const char* argv[])
{
	if (argc != 2) {
		fprintf(stderr, "Usage: %s input_dir", argv[1]);
		exit(1);
	}
}

int main(const int argc, const char* argv[])
{
	check_args(argc, argv);

	System* system   = newSystem();
	Parameter* param = newParameter(argv[1]);
	readParameterFromFile(param);
	Boundary* boundary = newBoundary(getBoundaryName(param));
	if (getBoundaryType(boundary) == PERIODIC) setBoxLength(boundary, getBoxlength(param));

#ifdef SIMULATION_3D
		initializeSystem(system, boundary, param, createFlatMesh, newTopolMesh);
#else
		if (getBoundaryType(boundary) == PERIODIC) {
			initializeSystem(system, boundary, param, createStraightChain, newTopolChain);
		} else if (getBoundaryType(boundary) == FREE) {
			initializeSystem(system, boundary, param, createRandomChain, newTopolChain);
		}
#endif

	// readRestartConfig(system, param);
	executeSimulation(system, boundary, param);
	// writeFinalConfig(system, param);

	deleteSystem(system);
	deleteParameter(param);
	deleteBoundary(boundary);
	return 0;
}
