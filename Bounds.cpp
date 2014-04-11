#include <string>
#include <iostream>
#include "Bounds.h"
using namespace std;

double Bounds::span(string dim) {
	if (dim == "x") {
		return xhi - xlo;
	} else if (dim == "y") {
		return yhi - ylo;
	} 
	cout << "\nBad Dim for bound span\n";
	return 0;
}

void Bounds::resize(string dim, double mult, double around) {
	double origin;
	if (dim =="x") {
		origin = xlo + around * (xhi - xlo);
		xhi = mult * (xhi - origin) + origin;
		xlo = mult * (xlo - origin) + origin;
	} else if (dim == "y") {
		origin = ylo + around * (yhi - ylo);
		yhi = mult * (yhi - origin) + origin;
		ylo = mult * (ylo - origin) + origin;
	}
}

bool Bounds::atomInBounds(Atom *a) {
	return a->pos.x >= xlo && a->pos.x <= xhi && a->pos.y >= ylo && a->pos.y <= yhi;
}
