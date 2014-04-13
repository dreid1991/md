#include "Python.h"
#include <stdio.h>
#include <iostream>
#include <string>
#include "Box.h"
#include "PythonTools.h"
#include "Atom.h"
#include "AtomGrid.h"
#include "AtomTools.h"
//#include "DataTools.h"
#include "list_macro.h"
#include "geometry_defs.h"
#include <regex>


using namespace std;
PythonTools pyTools;
DataTools dataTools;
AtomTools atomTools;

int main(int argc, const char *argv[]) {
	dataFileName = "lalala";
	InputData inputs = readInputData(dataFileName);
	double maxCutoff = inputs.maxCutoff;
	vector<Atom *> atoms = inputs.atoms;
	Box box = inputs.box;
	runNVE(atoms, box, 10, maxCutoff);
	return 0;
}
