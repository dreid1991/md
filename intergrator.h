#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <vector>
#include "math.h"
#include "Vector.h"
#include "AtomGrid.h"
void runInitStep(vector<Atom *> &, Box &);
void runStep(vector<Atom *> &, Box &);


#endif
