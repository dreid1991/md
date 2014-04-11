#ifndef BOUNDS_H
#define BOUNDS_H
#include "Vector.h"
#include <string>
#include "Atom.h"
using namespace std;
class Bounds {
	public:
		double xlo;
		double xhi;
		double ylo;
		double yhi;	
		Bounds() {
			xlo = 0;
			xhi = 0;
			ylo = 0;
			yhi = 0;
		}
		Bounds(double xlo_, double xhi_, double ylo_, double yhi_) {
			xlo = xlo_;
			xhi = xhi_;
			ylo = ylo_;
			yhi = yhi_;
		}
		double span(string dim);
		bool atomInBounds(Atom *);
		void resize(string dim, double mult, double around);

};

#endif
