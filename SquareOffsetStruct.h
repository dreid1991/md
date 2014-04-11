#ifndef SQUAREOFFSETSTRUCT_H
#define SQUAREOFFSETSTRUCT_H
#include <vector>
#include "Vector.h"
#include "Atom.h"
using namespace std;
struct square_offset {
	vector<vector<Atom *> *> gridSqrs;
	vector<Vector> offsets;
};

#endif
