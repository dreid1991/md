#include <iostream>
#include <vector>
#include "Box.h"
#include "Atom.h"
#include "Vector.h"
#include "PARMS.h"
#include "CreateAtoms.h"
using namespace std;

vector<Atom> CreateAtoms::box (Box &bounds, int count, int type) {
	map<int, SpeciesDef>::iterator defIt = PARMS::SPECIESDEFS.find(type);
	vector<Atom> atoms;
	if (defIt == SPECIESDEFS.end()) {
		cout << "BAD SPECIES DEF" << endl;
		return atoms;
	}
	SpeciesDef def = defIt->second;
	for (int i=0; i<count; i++) {
		atoms.push_back(Atom(def));
	}
	return atoms;
}


