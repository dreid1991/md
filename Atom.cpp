#include "Atom.h"
#include "Vector.h"


void Atom::applyForce(Vector &f) {
	a.x += f.x / m;
	a.y += f.y / m;
	a.z += f.z / m;
}

