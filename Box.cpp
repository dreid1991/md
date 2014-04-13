#include <string>
#include <iostream>
#include "Box.h"
using namespace std;


void Box::resize(string dim, double mult, double around) {
	double origin;
	if (dim =="x") {
		origin = xlo + around * (xhi - xlo);
		xhi = mult * (xhi - origin) + origin;
		xlo = mult * (xlo - origin) + origin;
	} else if (dim == "y") {
		origin = ylo + around * (yhi - ylo);
		yhi = mult * (yhi - origin) + origin;
		ylo = mult * (ylo - origin) + origin;
	} else if (dim == "z") {
		origin = zlo + around * (zhi - zlo);

		zhi = mult * (zhi - origin) + origin;
		zlo = mult * (zlo - origin) + origin;
	}
}

bool Box::atomInBox(Atom *a) {
	return a->pos.x >= xlo && a->pos.x <= xhi && a->pos.y >= ylo && a->pos.y <= yhi;
}
double Box::span(string dim) {
	if (dim == "x") {
		return xhi - xlo;
	} else if (dim == "y") {
		return yhi - ylo;
	} else if (dim == "z") {
		return zhi - zlo;
	}
	cout << "\nBad Dim for bound span\n";
	return 0;
}


void Box::loopIntoBox(Atom *a) {
	int nxLoop = floor((a->pos.x - xlo) / trace.x);
	int nyLoop = floor((a->pos.y - ylo) / trace.y); 
	int nzLoop = floor((a->pos.z - zlo) / trace.z);

	a->pos.x -= nxLoop * trace.x;
	a->pos.y -= nyLoop * trace.y;
	a->pos.z -= nzLoop * trace.z;
}
