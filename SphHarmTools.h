#ifndef SPHHARMTOOLS_H
#define SPHHARMTOOLS_H

#include <complex>
#include "Atom.h"
#include "Vector.h"
#include <vector>
#include "math.h"

class SphHarmTools {
	public:
		SphHarmTools() {

		}
		double PLM(int, int, double);
		int fact (int);
		complex<double> evalSphHarm(int, int, Vector &, int);
		vector<complex<double> > calcQlms(Atom *, int);
		double calcQl(Atom *, int);
		complex<double> calcQlmBar(Atom *, int, int);
		void calcQlAvg(vector<Atom *> &, int);
	
};

#endif
