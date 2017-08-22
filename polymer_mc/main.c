#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "system.h"
#include "topol.h"
#include "parameter.h"
#include "config_maker.h"
#include "boundary.h"

static void check_args(const int argc,
	const char* argv[])
{
	if (argc != 3) {
		fprintf(stderr, "Usage: %s input_dir [2D/3D]", argv[1]);
		exit(1);
	}
	const bool is3d = (strncmp(argv[2], "3D", 2) == 0);
	const bool is2d = (strncmp(argv[2], "2D", 2) == 0);
	if (!is3d && !is2d) {
		fprintf(stderr, "arg[2] should be [2D/3D].\n");
		exit(1);
	}
}

int main(const int argc, const char* argv[])
{
	check_args(argc, argv);
	const bool is3d = (strncmp(argv[2], "3D", 2) == 0);

	System* system   = newSystem();
	Parameter* param = newParameter(argv[1]);
	readParameterFromFile(param, is3d);
	Boundary* boundary = newBoundary(getBoundaryName(param));
	if (getBoundaryType(boundary) == PERIODIC) setBoxLength(boundary, getBoxlength(param));

	if (is3d) {
		initializeSystem(system, boundary, param, createFlatMesh, newTopolMesh);
	} else {
		initializeSystem(system, boundary, param, createStraightChain, newTopolChain);
	}

	// readRestartConfig(system, param);
	executeSimulation(system, boundary, param);
	// writeFinalConfig(system, param);

	deleteSystem(system);
	deleteParameter(param);
	deleteBoundary(boundary);
	return 0;
}