#include "runNVE.h"


void runNVE(vector<Atom *> &atoms, Box &box, int numTurns, double gridSize) {
	runInitStep(atoms, box);
	numTurns--;
	for (int i=0; i<numTurns; i++) {
		runStep(atoms, box);
	}
}


