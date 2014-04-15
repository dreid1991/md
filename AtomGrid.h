#ifndef ATOMGRID_H
#define ATOMGRID_H
#include "Atom.h"
#include <vector>
#include "Vector.h"
#include "Box.h"
#include "Grid.h"
#include "math.h"
#include <iostream>
#include <map>
#include "misc_math.h"

#include <string>
#include <algorithm>
#include "SquareOffsetStruct.h"
using namespace std;
class AtomGrid {
	public:
		Grid<vector<Atom *> > grid;
		Box box;
		double dx;
		double dy;
		double dz;
		unsigned int nx;
		unsigned int ny;
		unsigned int nz;
	//	vector<vector<Atom *> *> sqrs;
		vector<vector<vector<Atom *> > > &operator[](int x) {
			return grid[x];
		}
		Box getSquareBox(int, int int);

		map<map<string, double>, double> radii;
		Vector size;
		vector<Atom *> atoms;
		AtomGrid();	
		AtomGrid(vector<Atom *> &atoms, Box bounds, double dx, double dy, double dz, bool);
		AtomGrid(Box, double dx, double dy, double dz);
	//vector<vector<vector<Atom *> > > makeEmptyGrid(int nx, int ny);
		void printAtom(Atom *a);
		Vector sqrPosition(vector<Atom *> *);
		void populateAssignNeigh(vector<Atom *>, double cut);


};

#endif
