#include "Vector.h"
#include <vector>
#include "SpeciesDef.h"
#ifndef ATOM_H
#define ATOM_H

class Atom {
	public:
		Atom(SpeciesDef &);
		int type;
		int id;
		int sig;
		int eps;
		Vector p;
		Vector v;
		Vector a;
		double m;
		void applyForce(Vector &);
};
#endif
