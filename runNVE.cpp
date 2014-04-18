#include "runNVE.h"


void runNVE(vector<Atom *> &atoms, Box &box, int numTurns, double ljcut) {
	runInitStep(atoms, box, 1.2 * ljcut);
	numTurns--;
	for (int i=0; i<numTurns; i++) {
		runStep(atoms, box);
	}
}


