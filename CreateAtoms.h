#ifndef CREATEATOMS_H
#define CREATEATOMS_H
#include <vector>
#include "Atom.h"
#include "Box.h"

using namespace std;
class CreateAtoms {
	public:
		CreateAtoms();
		vector<Atom> box(Box &, int type, int count);
};

#endif
