#ifndef DATATOOLS_H
#define DATATOOLS_H
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include "math.h"
#include "Box.h"
#include "Atom.h"
#include <algorithm>
using namespace std;

string intToStr(int);

class DataTools {
	public:
		DataTools(){} 
		void binAvg(map<int, map<string, vector<double> > > *runData, string xPath, string yPath, vector<double> *xs, vector<double> *ys, int numBins, bool mvAvg, int span);
		double triangleAvg(vector<double> *xs, vector<double> *ys, int i, int span);
		vector<double> movingAvg(vector<double> *xs, vector<double> *ys, int span);
		double average(vector<double> *xs);
		vector<double> ddx(vector<double> *, vector<double> *, int order);
		void appendList(vector<double> *, vector<double> *);
		void hist(vector<double> *events, vector<double> *xs, vector<double> *ys, double xSpan, bool normalize);
		void histNoTri(vector<double> *events, vector<double> *xs, vector<double> *ys, double xSpan, bool normalize);
		void intHist(vector<int> *events, vector<double> *xs, vector<double> *ys, bool normalize);
		double linSlope(vector<double> &, vector<double> &);
};


#endif


