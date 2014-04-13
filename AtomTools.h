#ifndef ATOMTOOLS_H
#define ATOMTOOLS_H
#include "Python.h"
#include "Box.h"
#include "math.h"
#include "Atom.h"
#include <climits>
#include <vector>
#include <algorithm>
#include "Grid.h"
#include "PythonTools.h"
#include <iostream>
#include "Snap.h"
#include "misc_math.h"
#include <utility>
#include "displacement_datum.h"
#include "stamped_item.h"
#include "list_macro.h"
#include "DataTools.h"
#include "geometry_defs.h"
#include "timestamped_pos_def.h"
#include "MobilityProfile.h"
extern DataTools dataTools;
extern PythonTools pyTools;

class AtomTools {

	void calcMobilityProfile(MobilityProfile &, vector<Atom *> &, int);
	public:
		AtomTools(){}
		vector<Atom *> inBox(vector<Atom *> *atoms, Box b);
		vector<Atom *> inBox(vector<Atom> *atoms, Box b);
		void assignPositionHist(vector<Snap> &);
		Box findBounds(vector<Atom *> *);
		
		void assignNeighbors(vector<Atom *>, double, Box);
		void assignNeighsFromSquares(Grid<vector<Atom *> >, Atom *, int, int, int);
		vector<double> crystalGroupsTilt(vector<vector <Atom *> *> &crystalGroups);
		double crystalGroupTilt(vector<Atom *> &crystalGroup);
		
		vector<vector<int> > voronoiPolygonVertexIdxss(PyObject *voronoiData); 
		vector<int> voronoiTypeCount(vector<vector<int> > &);
		vector<int> voronoiTypeCount(vector<Atom * > &);
		vector<Polygon> voronoiPolygons(PyObject *, vector<vector<int> > &);
		
		vector<pair<double, displacement_datum> > calcDisplacementData(vector<Atom *> &);
		void avgStdevDisplacement(vector<pair<double, displacement_datum> > &disps); 
		vector<Atom *> selectByGrpSize(vector<vector<Atom *> *> &, int min, int max);
		vector<Atom *> test(const vector<Atom *> *, bool (*) (Atom *));
		int count(const vector<Atom *> *, bool (*) (Atom *));
		vector<vector<Atom *> *> assignCrystalGroups(vector<Atom *>, int, double);
		void growSeed(Atom *seed, vector<Atom *> *seeds, int qType, double qMin);
		void removeAtom(Atom *, vector<Atom *> *);
		double avgQ(vector<Atom *> *, int q);
		double testNeighborsFrac(vector<Atom *> *, bool(*) (Atom *), bool(*) (Atom *, Atom *));
		void assignAtomPositionHist(vector<Snap> &);
		void loopPositionHists(vector<Atom *> &, Box);
		void loopAtomHist(vector<timestamped_pos> &, unsigned int curAtomIdx, Box &, bool, bool);

		
		vector<stamped_item<double, MobilityProfile> > calcMobilityProfiles(vector<Atom *> &, double, double, double); 
};

#endif
