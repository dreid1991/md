#ifndef RUNNVE_H
#define RUNNVE_H

#include <vector>
#include "Box.h"
#include "Atom.h"
#include "Vector.h"
#include "AtomGrid.h"
using namespace std;

void runNVE(vector<Atom *> &, Box &, int);

#endif
