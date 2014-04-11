#ifndef MISCVDEP_H
#define MISCVDEP_H
#include "Python.h"
#include "AtomTools.h"
#include "experiment.h"
#include "DataTools.h"
#include "PythonTools.h"
#include "SphHarmTools.h"
#include "math.h"
#include <functional>
#include <numeric>
#include <stdio.h>
#include <vector>
#include <string>
#include <map>
#define CALCEXPSTDEV(e, src, dest)\
	calcExpStdev(e,\
	[&] (run& r) {return r.src;},\
	[&] (experiment& e, double val) {e.dest = val;}) 

#define MAKEATOMGROUP(leExp, leRun, testFunc, attrList, attrFrac)\
	makeAtomGroup(leExp, leRun, testFunc,\
	[&] (run& r, vector<Atom *>& list) {r.attrList = list;},\
	[&] (experiment& e, vector<Atom *>& list) {e.attrList.insert(e.attrList.end(), list.begin(), list.end());},\
	[&] (run& r, double frac) {r.attrFrac = frac;}\
	)

#define SETAGVAL(exp, fieldToSet, runVarName, fracExpr)\
	setAgVal(exp, \
	[&] (run& runVarName) {return fracExpr;},\
	[&] (experiment& e, double frac) {e.fieldToSet = frac;},\
	[&] (run& r, double frac) {r.fieldToSet = frac;}\
	)



extern double xrylThresh;
extern double amphThresh;
extern PythonTools pyTools;
extern DataTools dataTools;
extern AtomTools atomTools;
extern SphHarmTools sphTools;	

bool atomType(Atom *a);

bool type1(Atom *a);
bool type2(Atom *a);

double roundDouble(double, int);

bool alwaysTrue(Atom *a);

bool isXryl(Atom *a);

bool isAmph(Atom *a);

bool isTrans(Atom *a);
void setAgVal(experiment& e, function<double (run&)> getFrac, function<void (experiment&, double)> setExpVal, function<void (run&, double)> setRunVal);

void calcExpStdev(experiment &e, function<double (run&)> fn, function<void (experiment&, double)> setVal); 

bool diffType(Atom *a, Atom *n);

void eraseExt(map<int, vector<string> > &fnsMap);

void assignQGroups(map<int, experiment> &exps);

void populateQ6s(map<int, experiment> &exps);

void makeTitles(map<int, experiment> &exps); 

void processGrps(map<int, vector<string> > &fnsBase, map<int, experiment> &exps);

void assignXrylGroups(map<int, experiment> &exps); 


void convertToPys(map<int, experiment> &exps, PyObject *xssPy, PyObject *yssPy, PyObject *titlesPy);


void convertToPysXrylSize(map<int, experiment> &exps, PyObject *xss, PyObject *yss, PyObject *titlesPy);

void assignXrylHist(map<int, experiment> &exps); 

void push_backExps(map<int, experiment> &exps, vector<experiment *> &all);

#endif
