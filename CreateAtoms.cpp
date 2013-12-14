#include <iostream>
#include <vector>
#include "Box.h"
#include "Atom.h"
#include "Vector.h"
#include "PARMS.h"
#include "CreateAtoms.h"
#include "SpeciesDef.h"

#include <stdlib.h>
using namespace std;
CreateAtoms::CreateAtoms() {

}
vector<Atom> CreateAtoms::box (Box &bounds, int type, int count) {
	SpeciesDef foo = PARMS::SPECIESDEFS.find(1)->second;
	cout << foo.m << endl;
	map<int, SpeciesDef>::iterator defIt = PARMS::SPECIESDEFS.find(type);
	vector<Atom> atoms;

	if (defIt == PARMS::SPECIESDEFS.end()) {
		cout << "BAD SPECIES DEF" << endl;
		return atoms;
	}
	SpeciesDef def = defIt->second;
	for (int i=0; i<count; i++) {
		Vector p = bounds.randV();
		cout << p.x << endl	;
		atoms.push_back(Atom(def, p));
	}
	return atoms;
}


