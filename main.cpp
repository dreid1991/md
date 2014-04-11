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
	vector<Atom *> atoms = inputs.atoms;
	Box box = inputs.box;
	return 0;
}
