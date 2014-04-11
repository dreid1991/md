#ifndef DISPLACEMENT_DATUM_H
#define DISPLACEMENT_DATUM_H

#include <map>
#include <vector>

using namespace std;

struct displacement_datum {
	vector<double> displacements;
	double avg;
	double stderror;
	double stdev; 
};

displacement_datum displacement_datum_init();

#endif
