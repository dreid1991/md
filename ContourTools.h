#ifndef CONTOURTOOLS_H
#define CONTOURTOOLS_H

#include "Python.h"
#include "Atom.h"
#include <vector>
#include "PythonTools.h"
#include <map>
#include "Vector.h"
#include "list_macro.h"
#include "SquareOffsetStruct.h"
#include "AtomGrid.h"
#include "misc_math.h"
class ContourTools {
	PythonTools *pyTools;
	public:
		ContourTools(PythonTools *pyTools_) {
			pyTools = pyTools_;

		}
		void graphQContour(vector<Atom *> &atoms, AtomGrid &grid, int qType, double resolution, bool loopX, bool loopyY, double figSize, double searchRad, map<string, double>, int numLevels, double maxLevel, bool scatterAtoms);
		bool pointIsAlone(double x, double y, square_offset *square_offsets, double searchRad, map<string, double> relRadii);
		void assignZ(square_offset *square_offsets, double x, double y, double *zPtr, int qType); 
};

#endif

