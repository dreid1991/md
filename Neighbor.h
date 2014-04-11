#ifndef NEIGHBOR_H
#define NEIGHBOR_H
#include "Atom.h"
#include "Vector.h"
class Neighbor {
	public:
		Neighbor(Atom *a, double offX, double offY): offset(Vector(offX, offY)) {

		}
		Atom *atom;
		Vector offset;
};



#endif 
