#ifndef BOUNDS_H
#define BOUNDS_H
#include "Vector.h"
#include <string>
#include "Atom.h"
using namespace std;
class Box {
	public:
		double xlo;
		double xhi;
		double ylo;
		double yhi;	
		double zlo;
		double zhi;
		Box() {
			xlo = 0;
			xhi = 0;
			ylo = 0;
			yhi = 0;
			zlo = 0;
			zhi = 0;
		}
		Box(double xlo_, double xhi_, double ylo_, double yhi_, double zlo_, double zhi_) {
			xlo = xlo_;
			xhi = xhi_;
			ylo = ylo_;
			yhi = yhi_;
			zlo = zlo_;
			zhi = zhi_;
			trace = Vector(xhi - xlo, yhi - ylo, zhi - zlo);
		}
		Vector trace;
		double span(string dim);
		bool atomInBox(Atom *);
		void resize(string dim, double mult, double around);
		void loopIntoBox(Atom *);

};

#endif
