#ifndef EXPERIMENT_H
#define EXPERIMENT_H
#include <vector>
#include <map>
#include <string>
#include "Atom.h"
#include "AtomGrid.h"
#include "displacement_datum.h"
using namespace std;

struct run {
	vector<Atom *> atoms;
	vector<Atom *> atomsInBounds;
	vector<Atom *> atomsXryl;
	vector<Atom *> atomsAmph;
	vector<Atom *> atomsTrans;
	vector<Atom *> atomsAmphTrans;
	double frac1Bulk;
	double frac1Xryl;
	double frac1Amph;
	double frac1Trans;
	double frac1AmphTrans;

	double fracDiffTypeXryl;
	double fracDiffTypeAmph;
	double fracDiffTypeBulk;	
	double fracDiffTypeAmphTrans;	
	
	double fracXryl;
	double fracAmph;
	double fracTrans;
	double fracAmphTrans;

	double ratioType1XrylAmph;

	AtomGrid grid;
	double eng;
	vector<vector<Atom *> *> crystalGroups;
};
struct experiment {
	vector<run> runs;	
	vector<double> xssQ6Hist; 
	vector<double> yssQ6Hist;
	vector<double> xsXrylSizeHist;
	vector<double> ysXrylSizeHist;

	map<double, displacement_datum> displacementData; 	

	vector<vector<Atom *> *> crystalGroups;
	double frac1LrgGrps;
	double frac1LrgGrpsStdev;
	double frac1SmlGrps;
	string title;
	double qXrylAvg;
	double qAmphAvg;
	double qTransAvg;

	double fracXryl;
	double fracAmph;
	double fracTrans;
	double fracAmphTrans;

	double fracXrylStdev;
	double fracAmphStdev;
	double fracTransStdev;
	double fracAmphTransStdev;

	double frac1Xryl;
	double frac1Amph;
	double frac1Bulk;
	double frac1AmphTrans;

	double ratioType1XrylAmph;
	double ratioType1XrylAmphStdev;

	double frac1XrylStdev;
	double frac1AmphStdev;
	double frac1BulkStdev;
	double frac1AmphTransStdev;

	double fracDiffTypeXryl;
	double fracDiffTypeAmph;
	double fracDiffTypeBulk;
	double fracDiffTypeAmphTrans;


	double fracDiffTypeXrylStdev;
	double fracDiffTypeAmphStdev;
	double fracDiffTypeBulkStdev;
	double fracDiffTypeAmphTransStdev;

	int id;
	string idStr;
	vector<Atom *> atomsXryl;
	vector<Atom *> atomsAmph;
	vector<Atom *> atomsTrans;
	vector<Atom *> atomsAmphTrans;
	vector<Atom *> atomsInBounds;
	vector<Atom *> type12;
};
#endif

