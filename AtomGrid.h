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
		square_offset getNeighborSquares(int initX, int initY, bool loopX, bool loopY);
		Box box;
		double dx;
		double dy;
		unsigned int nx;
		unsigned int ny;
		
		vector<vector<Atom *> > &operator[](int x) {
			return grid[x];
		}
		double getMaxY(vector<Atom *> &);
		double avgSurfaceY();

		map<map<string, double>, double> radii;
		void sqrIdx(int *, int *, double x, double y); 
		Vector size;
		vector<Atom *> *atoms;
		double findRadius(map<string, double> relRadii);
		AtomGrid();	
		AtomGrid(vector<Atom *> &atoms, Box bounds, double dx, double dy);
	//vector<vector<vector<Atom *> > > makeEmptyGrid(int nx, int ny);
		void assignNeighborsFromSqr(Atom *a, vector<Atom*>*sqr, Vector offset, double rSqr);
		void assignNeighbors(double rThresh, bool loopX, bool loopY);
		void assignNeighborsAtom(Atom *a, square_offset *sqrOffs, double rSqr);
		void printAtom(Atom *a);
		Vector sqrPosition(vector<Atom *> *);


};

#endif
