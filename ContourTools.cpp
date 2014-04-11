#include "ContourTools.h"
#include "math.h"
#include "PythonTools.h"
#include <iostream>
void ContourTools::assignZ(square_offset *square_offsets, double x, double y, double *zPtr, int qType) {
	double weight = 0;
	Vector pos = Vector(x, y);
	for (unsigned int i=0; i<square_offsets->gridSqrs.size(); i++) {
		vector<Atom *> *gridSqr = square_offsets->gridSqrs[i];
		Vector offset = square_offsets->offsets[i];
		for (unsigned int nIdx=0; nIdx<gridSqr->size(); nIdx++) {
			double dist = fmax(1e-5, pos.dist((*gridSqr)[nIdx]->pos + offset));

			(*zPtr) += (*gridSqr)[nIdx]->qs[qType] / (dist * dist);
			weight += 1 / (dist * dist);
		}
	}
	(*zPtr) /= weight;
}

bool ContourTools::pointIsAlone(double x, double y, square_offset *square_offsets, double searchRad, map<string, double> relRadii) {
	Vector gridPt = Vector(x, y);
	double radSqr = searchRad * searchRad;
	for (unsigned int i=0; i<square_offsets->gridSqrs.size(); i++) {
		vector<Atom *> *gridSqr = square_offsets->gridSqrs[i];
		Vector *offset = &square_offsets->offsets[i];
		for (unsigned int nIdx=0; nIdx<gridSqr->size(); nIdx++) {
			Atom *a = (*gridSqr)[nIdx];
			Vector aPt = Vector(a->pos.x + offset->x, a->pos.y + offset->y);
			if (aPt.distSqr(gridPt) < radSqr * relRadii[a->type] * relRadii[a->type]) {
				return false;
			}
		}
	}
	return true;
}


void ContourTools::graphQContour(vector<Atom *> &atoms, AtomGrid &grid, int qType,  double resolution, bool loopX, bool loopY, double figSize, double searchRad, map<string, double> relRadii, int numLevels, double maxLevel, bool scatterAtoms) {
	if (searchRad == -1) {
		searchRad = 1.5 * grid.findRadius(relRadii); 
	}
	double maxQVal = 0;
	for (unsigned int i=0; i<atoms.size(); i++) {
		maxQVal = fmax(atoms[i]->qs[qType], maxQVal);
	}
	int numX = (int) ceil(grid.bounds.span("x") / resolution) + 1;
	int numY = (int) ceil(grid.bounds.span("y") / resolution) + 1;
	vector<double> xs; xs.reserve(numX);
	vector<double> ys; ys.reserve(numY);
	for (int i=0; i<numX; i++) {
		xs.push_back(i * resolution);
	}
	for (int i=0; i<numY; i++) {
		ys.push_back(i * resolution);
	}
	Grid<double> zs = Grid<double>(numY, numX);//row, col
	for (int xIdx=0; xIdx<numX; xIdx++) {
		double x = grid.bounds.xlo + xIdx * resolution;
		for (int yIdx=0; yIdx<numY; yIdx++) {

			double y = grid.bounds.ylo + yIdx * resolution;
			int xGridIdx;
			int yGridIdx;
			grid.sqrIdx(&xGridIdx, &yGridIdx, x, y);

			square_offset square_offsets = grid.getNeighborSquares(xGridIdx, yGridIdx, loopX, loopY); 
			if (!pointIsAlone(x, y, &square_offsets, searchRad, relRadii)) {
				assignZ(&square_offsets, x, y, &zs[yIdx][xIdx], qType);
			}
		}
	}
	vector<double> levels;
	if (numLevels != -1 && maxLevel != -1) {
		for (int i=0; i<numLevels + 1; i++) {
			double level = maxLevel * i / numLevels;
			levels.push_back(roundDouble(level, 3));
		}
	}
	double xSize = grid.bounds.span("x");
	double ySize = grid.bounds.span("y");
	double ratio = xSize / figSize;
	xSize /= ratio;
	ySize /= ratio;
	xSize *= 1.2;
	PyObject *graphMod = PyImport_ImportModule("graphQContour");
	PyObject *graphFunc = PyObject_GetAttrString(graphMod, "graphQContour");
	
	PyObject *atomXs;
	PyObject *atomYs;
	if (scatterAtoms) {
		atomXs = LISTMAPPY(Atom *, a, atoms, PyFloat_FromDouble(a->pos.x));
		atomYs = LISTMAPPY(Atom *, a, atoms, PyFloat_FromDouble(a->pos.y));
	}


	PyObject *args = PyTuple_New(12);
	PyObject *xFloat = PyFloat_FromDouble(xSize);
	PyObject *yFloat = PyFloat_FromDouble(ySize);

	PyObject *pyZs = PyList_New(zs.size());

	for (unsigned int i=0; i<zs.size(); i++) {
		PyObject *listAsPy = pyTools->doubleVecToPy(&zs[i]);
		PyList_SetItem(pyZs, i, listAsPy);
	//	Py_DECREF(listAsPy);
	}

	PyObject *pyXs = pyTools->doubleVecToPy(&xs);
	PyObject *pyYs = pyTools->doubleVecToPy(&ys);
	PyObject *pyLevels = pyTools->doubleVecToPy(&levels);

	
	char buffer[5];
	sprintf(buffer, "6%d", qType);
	PyObject *label = PyString_FromString(buffer);
	PyTuple_SetItem(args, 0, pyTools->plt);
	
	PyTuple_SetItem(args, 1, pyXs);
	PyTuple_SetItem(args, 2, pyYs);
	PyTuple_SetItem(args, 3, pyZs);
	PyTuple_SetItem(args, 4, pyLevels);
	PyTuple_SetItem(args, 5, xFloat);
	PyTuple_SetItem(args, 6, yFloat);
	PyTuple_SetItem(args, 7, PyFloat_FromDouble(maxLevel));
	PyTuple_SetItem(args, 8, PyFloat_FromDouble(maxQVal));
	PyTuple_SetItem(args, 9, label);
	PyTuple_SetItem(args, 10, atomXs);
	PyTuple_SetItem(args, 11, atomYs);
	PyObject_CallObject(graphFunc, args);
	Py_DECREF(pyXs);
	Py_DECREF(pyYs);
	Py_DECREF(pyZs);
	Py_DECREF(xFloat);
	Py_DECREF(yFloat);
	Py_DECREF(pyLevels);
	Py_DECREF(graphMod);
	Py_DECREF(graphFunc);
	Py_DECREF(atomXs);
	Py_DECREF(atomYs);
	Py_DECREF(args);

}

