#include "Vector.h"
#include <vector>

#ifndef ATOM_H
#define ATOM_H

class Atom {
	public:
		int type;
		int id;
		Vector p;
		Vector v;
		Vector a;
		double m;
		void applyForce(Vector &);
};
#endif
