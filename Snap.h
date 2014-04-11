#ifndef SNAP_H
#define SNAP_H
#include "Atom.h"
#include <vector>
#include "Bounds.h"
using namespace std;
struct Snap {
	vector<Atom *> atoms;
	Bounds bounds;
	double timestamp;
};
#endif
