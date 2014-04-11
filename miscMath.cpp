#include "misc_math.h"
using namespace std;
double roundDouble(double x, int numDec) {
	double y = pow(10, numDec);
	x *= y;
	x = (double) floor(x + .5);
	x /= y;
	return x;
}

double average(vector<double> &xs) {
	double sum = 0;
	for (unsigned int i=0, ii=xs.size(); i<ii; i++) {
		sum += xs[i];
	}
	return sum / xs.size();
}


