#include "Atom.h"
#include "Vector.h"

Atom::Atom(SpeciesDef &def, Vector p_): p(Vector(p_)), v(Vector()), a(Vector()) {
	m = def.m;
	type = def.type;
	sig = def.sig;
	eps = def.eps;
}
void Atom::applyForce(Vector &f) {
	a.x += f.x / m;
	a.y += f.y / m;
	a.z += f.z / m;
}

