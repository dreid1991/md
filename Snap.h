#ifndef SNAP_H
#define SNAP_H
#include "Atom.h"
#include <vector>
#include "Box.h"
using namespace std;
struct Snap {
	vector<Atom *> atoms;
	Box box;
	double timestamp;
};
#endif
