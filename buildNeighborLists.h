#ifndef BUILDNEIGHBORLISTS_H
#define BUILDNEIGHBORLISTS_H
#include "math.h"
#include "AtomGrid.h"
#include "Atom.h"

#include "Box.h"
#include <vector>
void buildNeighborLists(vector<Atom *> &, Box, int, int, int, double ljCut, double neighCut);

#endif
