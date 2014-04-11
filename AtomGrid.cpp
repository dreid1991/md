
#include "AtomGrid.h"
using namespace std;


AtomGrid::AtomGrid(vector<Atom *> *atoms_, Bounds bounds_, double dx_, double dy_): bounds(bounds_) {
	atoms = atoms_;
	dx = dx_;
	dy = dy_;
	nx = (unsigned int) ceil((bounds.xhi - bounds.xlo) / dx);
	ny = (unsigned int) ceil((bounds.yhi - bounds.ylo) / dy);
	size = Vector(bounds.xhi - bounds.xlo, bounds.yhi, bounds.ylo);
	grid = Grid<vector<Atom *> > (nx, ny);
	for (unsigned int i=0; i<atoms->size(); i++) {
		unsigned int x = (int) (*atoms)[i]->pos.x / dx;
		unsigned int y = (int) (*atoms)[i]->pos.y / dy;
		if (x < grid.size() && y < grid[0].size()) {
			grid[x][y].push_back((*atoms)[i]);
		} else if ((*atoms)[i]->pos.x == size.x || (*atoms)[i]->pos.y == size.y) {
			cout << "Looping atom from " << (*atoms)[i]->pos.x << ", " << (*atoms)[i]->pos.y << endl;
			(*atoms)[i]->pos.x -= size.x * (int) ((*atoms)[i]->pos.x / size.x);
			(*atoms)[i]->pos.y -= size.y * (int) ((*atoms)[i]->pos.y / size.y);
			cout << "Looped atom to " << (*atoms)[i]->pos.x << ", " << (*atoms)[i]->pos.y << endl;
		} else {
			cout << "Atom off grid" << endl;
			(*atoms)[i]->print();
		}
	}
}

AtomGrid::AtomGrid() {

}

double AtomGrid::getMaxY(vector<Atom *> &as) {
	double maxY = -10000;
	for (unsigned i=0; i < as.size(); i++) {
		maxY = max(maxY, as[i]->pos.y);
	}
	return maxY;
}


double AtomGrid::avgSurfaceY() {
	vector<double> maxYs;
	for (int x=0; x<grid.nx; x++) {
		for (int y = grid.ny-1; y>=0; y--) {
			if (grid[x][y].size()) {
				maxYs.push_back(getMaxY(grid[x][y]));
				cout << maxYs[maxYs.size() - 1] << endl;
				break;
			}
		}
	}
	return average(maxYs);
}

Vector AtomGrid::sqrPosition(vector <Atom*> *sqr) {
	for (unsigned int i=0; i<grid.size(); i++) {
		for (unsigned int j=0; j<grid[i].size(); j++) {
			if (&grid[i][j] == sqr) {
				return Vector(bounds.xlo + dx * i, bounds.ylo + dy * j);
			}
		}
	}
	return Vector(-1, -1);
}

void AtomGrid::printAtom(Atom *a) {
	cout << "\nPrintingAtom\n";
	cout << "id: ";
	cout << a->id;
	cout << "\ntype: ";
	cout << a->type;
	cout << "\nx: ";
	cout << a->pos.x;
	cout << "\ny: ";
	cout << a->pos.y;
	cout << "\n";
}



double AtomGrid::findRadius(map<string, double> relRadii) {
	if (radii.find(relRadii) == radii.end()) {
		const unsigned int step = 2; //because who needs to check all of them?
		double radius = 0;
		int numAtomsUsed = 0;
		for (unsigned int i=0; i<atoms->size(); i+=step) {
			Atom *a = (*atoms)[i];
			vector<Atom *> *neighbors = &a->neighbors;
			vector<Vector> *offsets = &a->neighborOffsets;
			double minNeighborRad = INFINITY;
	
			for (unsigned int nIdx=0; nIdx<neighbors->size(); nIdx++) {
				double dist = a->pos.dist((*neighbors)[nIdx]->pos + (*offsets)[nIdx]);
				double rad = dist / (relRadii[a->type] + relRadii[(*neighbors)[nIdx]->type]);
				minNeighborRad = fmin(minNeighborRad, rad);
			}
			if (minNeighborRad != INFINITY) {
				radius = radius * numAtomsUsed / (numAtomsUsed + 1) + minNeighborRad / (numAtomsUsed + 1);
			}
		}
		radii[relRadii] = radius;
		return radius;
	} else {
		return radii[relRadii];
	}
}

square_offset AtomGrid::getNeighborSquares(int initX, int initY, bool loopX, bool loopY) {
	vector<vector<Atom *> *> gridSqrs;
	vector<Vector> offsets;
	for (int x=initX-1; x<=initX+1; x++) {
		for (int y=initY-1; y<=initY+1; y++) {
			int flrX = 0;
			int flrY = 0;
			if (x < 0) {
				flrX = -1;
			} else if (x >= (signed int) nx) {
				flrX = 1;
			}
			if (y < 0) {
				flrY = -1;
			} else if (y >= (signed int) ny) {
				flrY = -1;
			}

			bool appendMe;
			if (flrX && flrY) {
				appendMe = loopX && loopY;
			} else if (flrY) {
				appendMe = loopY;
			} else if (flrX) {
				appendMe = loopX;
			} else {
				appendMe = true;
			}
			if (appendMe) {
				double offsetX = flrX * size.x;
				double offsetY = flrY * size.y;
				int offsetIdxX = x - flrX * nx;
				int offsetIdxY = y - flrY * ny;
				gridSqrs.push_back(&grid[offsetIdxX][offsetIdxY]);
				Vector offset = Vector(offsetX, offsetY);

				offsets.push_back(offset);
			}
		}
	}
	square_offset toReturn;
	toReturn.gridSqrs = gridSqrs;
	toReturn.offsets = offsets;
	return toReturn;
}

void AtomGrid::assignNeighborsFromSqr(Atom *a, vector<Atom *> *gridSqr, Vector offset, double rSqr) {
	int len = gridSqr->size();
	for (int i=0; i<len; i++) {
		Atom *neighbor = (*gridSqr)[i];
		Vector neighborPos = Vector(neighbor->pos.x + offset.x, neighbor->pos.y + offset.y);
		if (a->pos.distSqr(neighborPos) < rSqr && a != neighbor) {
			a->neighbors.push_back(neighbor);
			a->neighborOffsets.push_back(offset);
		}
	}
}

void AtomGrid::sqrIdx(int *xIdx, int *yIdx, double x, double y) {
	(*xIdx) = (int) (((x - bounds.xlo) / dx));
	(*yIdx) = (int) (((y - bounds.ylo) / dy));
}

void AtomGrid::assignNeighborsAtom(Atom *a, square_offset *sqrOffs, double rSqr) {
	vector<vector<Atom *> *> *gridSqrs = &(*sqrOffs).gridSqrs;
	vector<Vector> *offsets = &(*sqrOffs).offsets;
	for (unsigned int i=0; i<gridSqrs->size(); i++) {
		assignNeighborsFromSqr(a, (*gridSqrs)[i], (*offsets)[i], rSqr);
	}
}

void AtomGrid::assignNeighbors(double rThresh, bool loopX, bool loopY) {
	double rSqr = rThresh * rThresh;
	for (unsigned int x=0; x<nx; x++) {
		for (unsigned int y=0; y<ny; y++) {
			vector<Atom *> *gridSqr = &grid[x][y];
			square_offset sqrOffs = getNeighborSquares(x, y, loopX, loopY);
			for (unsigned int i=0; i<gridSqr->size(); i++) {
				assignNeighborsAtom((*gridSqr)[i], &sqrOffs, rSqr);
			}
		}
	}
}



