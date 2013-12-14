#ifndef CREATEATOMS_H
#define CREATEATOMS_H
#include <vector>
#include "Atom.h"
#include "Box.h"

using namespace std;
class CreateAtoms {
	public:
		vector<Atom> box(Box &, int count, int type);
};

#endif
