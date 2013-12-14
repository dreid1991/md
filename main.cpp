#include <iostream>
#include <vector>
#include "Vector.h"
#include "Atom.h"
#include "Box.h"
#include <map>
#include "CreateAtoms.h"
#include "SPECIESDEFS.h"
#include "SpeciesDef.h"
#include "PARMS.h"
using namespace std;
map<int, SpeciesDef> PARMS::SPECIESDEFS;

int main() {
	asserts::assertSpeciesDefs();
	CreateAtoms createAtoms = CreateAtoms();
	Vector origin = Vector(0, 0, 0);
	Vector trace = Vector(10, 10, 10);
	Box bounds = Box(origin, trace);
	vector<Atom> atoms = createAtoms.box(bounds, 1, 10);
	for (unsigned int i=0; i<atoms.size(); i++) {
		cout << atoms[i].m << endl;
	}
}
