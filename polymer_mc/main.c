#include <stdio.h>
#include <stdlib.h>

#include "system.h"
#include "topol.h"
#include "parameter.h"
#include "config_maker.h"

int main(int argc, char* argv[])
{
	if (argc != 2) {
		fprintf(stderr, "Usage: %s input_dir", argv[1]);
		exit(1);
	}
	System* system   = newSystem();
	Parameter* param = newParameter(argv[1]);
	readParameterFromFile(param);
	initializeSystem(system, param, createStraightChain, newTopolChain);

	// readRestartConfig(system, param);
	executeSimulation(system, param);
	// writeFinalConfig(system, param);

	deleteSystem(system);
	deleteParameter(param);
	return 0;
}