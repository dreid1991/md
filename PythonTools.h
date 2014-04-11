#ifndef PYTHONTOOLS_H
#define PYTHONTOOLS_H
#include "Python.h"
#include "Snap.h"
#include "DataTools.h"
#include "math.h"
#include <algorithm>
#include <stdlib.h>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sys/dir.h>
//#include <regex>
#include <dirent.h>
#include <sys/types.h>
#include "fileSpec.h"
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <error.h>
#include "list_macro.h"
using namespace std;
class PythonTools {
	public:
		PyObject *dumpFunc;
		PyObject *plt;
		PyObject *re;
		PyObject *reCompile;
		PyObject *reSearch;
		PyObject *isNone;
		PyObject *show;
		PyObject *savefig;
		PyObject *scatter;
		PyObject *plot;
		PyObject *grabDataMod;
		PyObject *hist;
		
		//scipy stuff
		PyObject *scipy;
		PyObject *spatial;
		PyObject *voronoi;
		PythonTools(){};
		
		//end scipy
		//Dump grabDump(string fn);
		string vectorToStr(vector<double> *);
		
		PyObject *callVoronoi(vector<Atom *> &);
		
		PyObject *convertNPArray(PyObject *, bool);

		void init();
		void setDump(PyObject *);
		void setPlt(PyObject *);
		void setGrab(PyObject *);
		bool allNums(string str);
		int countDataCols(vector<string> *);
		PyObject * doubleVecToPy(vector<double> *);
		void appendSnap(PyObject *pySnps, vector<Snap> *cSnps, int idx);
		vector<Snap> loadDump(string fn, double timestep);
		map<int, string> readDescData(string folderName, string entry);


		map<int, PyObject *> compileRes(map<int, string> &);
		map<int, vector<string> > grabFileNames(string folderName, map<int, string>);
		PyObject *mapIntDoubleToPy(map<int, double> &);
		map<string, vector<double> > grabData(string, long);
		vector<string> PyToVecString(PyObject *);
		vector<double> PyToVecDouble(PyObject *);
		void scatterAtoms(vector<Atom *> &, PyObject *dict);
		vector<string> grabVDepByDesc(vector<string> folders, string descItem, string val);

		vector<map<string, double> > grabOutComments(string folderName, string fileName, vector<string> comments);
};
#endif
