
#include "AtomGrid.h"
using namespace std;


AtomGrid::AtomGrid(vector<Atom *> &atoms_, Box box_, double dx_, double dy_, double dz_, bool atomsAreOnGrid): box(box_) {
	atoms = atoms_;
	dx = dx_;
	dy = dy_;
	dz = dz_;
	nx = (unsigned int) ceil((box.xhi - box.xlo) / dx);
	ny = (unsigned int) ceil((box.yhi - box.ylo) / dy);
	nz = (unsigned int) ceil((box.zhi - box.zlo) / dz);
	size = Vector(box.xhi - box.xlo, box.yhi - box.ylo, box.zhi - box.zlo);
	grid = Grid<vector<Atom *> > (nx, ny, nz);
/*	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			for (unsigned int k=0; k<nz; k++) {
				sqrs.push_back(&grid[i][j][k]);
			}
		}
	}*/
	//plan: make dense grid (one with 2 * cutoff spacing)
	//make grid with spacing of CPU sectors.  dense grid cells must fall exactly into one CPU grid sector
	//populate CPU sector grid with atoms.  
	//Take these CPU sector atom lists and populate the dense grid.  Add one atom at a time and add to neighbor lists.  
	//Do not deal with ghost particles in this step
	//Let all finish gridding
	//Now we must build ghost lists
	//for now, have each CPU run along border cells and add neighbors to its own atoms.
	if (atomsAreOnGrid) {
		for (unsigned int i=0; i<atoms.size(); i++) {
			Atom *a = atoms[i];

			int x = (a->pos.x - box.xlo) / box.trace.x;
			int y = (a->pos.y - box.ylo) / box.trace.y;
			int z = (a->pos.z - box.zlo) / box.trace.z;
			grid[x][y][z].push_back(a);
		}
	} else {
		for (unsigned int i=0; i<atoms.size(); i++) {
			Atom *a = atoms[i];
			int x = floor((a->pos.x - box.xlo) / box.trace.x);
			int y = floor((a->pos.y - box.ylo) / box.trace.y);
			int z = floor((a->pos.z - box.zlo) / box.trace.z);
			int dx = -floor(x / nx);
			int dy = -floor(y / ny);
			int dz = -floor(z / nz);
			a->pos.x += dx * box.trace.x;
			a->pos.y += dy * box.trace.y;
			a->pos.z += dz * box.trace.z;
			grid[x][y][z].push_back(a);
		}
	}
}
AtomGrid::AtomGrid() {

}
AtomGrid::AtomGrid(Box box, int dx_, int dy_, int dz_) {
	dx = dx_;
	dy = dy_;
	dz = dz_;
	nx = (unsigned int) ceil((box.xhi - box.xlo) / dx);
	ny = (unsigned int) ceil((box.yhi - box.ylo) / dy);
	nz = (unsigned int) ceil((box.zhi - box.zlo) / dz);
	size = Vector(box.xhi - box.xlo, box.yhi - box.ylo, box.zhi - box.zlo);
	grid = Grid<vector<Atom *> > (nx, ny, nz);
/*	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			for (unsigned int k=0; k<nz; k++) {
				sqrs.push_back(&grid[i][j][k]);
			}
		}
	}*/
}

void AtomGrid::populateAssignNeigh(vector<Atom *> &atoms, double cut) {
	double cutSqr = cut * cut;
	for (unsigned int i=0, ii=atoms.size(); i<ii; i++) {
		Atom *a = atoms[i];
		a->neighbors.erase(a->neighbors.begin(), a->neighbors.end()); //keeping memory alloced, why not.  who knows what's faster
		Vector aPos = a->pos;
		a->posWhenGrid = aPos;
		int x = (a->pos.x - box.xlo) / dx;
		int y = (a->pos.y - box.ylo) / dy;
		int z = (a->pos.z - box.zlo) / dz;
		for (int i = fmax(0, x-1), ii = fmin(nx, x+2); i<ii; i++) {
			for (int j = fmax(0, y-1), jj = fmin(ny, y+2); j<jj; j++) {
				for (int k = fmax(0, z-1), kk = fmin(nz, z+2); k<kk; k++) {
					vector<Atom *> sqr = grid[i][j][k];
					for (unsigned int q=0; q<sqr.size(); q++) {
						Atom *b = sqr[i];
						if (b->pos.distSqr(a) <= cutSqr) {
							a->neighbors.push_back(b);
							b->neighbors.push_back(a);
						}
					}
				}
			}
		}
	}
}

Vector AtomGrid::sqrPosition(vector <Atom*> *sqr) {
	for (unsigned int i=0; i<nx; i++) {
		for (unsigned int j=0; j<ny; j++) {
			for (unsigned int k=0; k<nz; k++) {
				if (&grid[i][j][k] == sqr) {
					return Vector(box.xlo + dx * i, box.ylo + dy * j, box.zlo + dz * k);
				}
			}
		}
	}
	return Vector(-1, -1, -1);
}

Box getSquareBox(int x, int y, int z) {
	return Box(box.xlo + x * dx, box.ylo + y * dy, box.zlo + z * dz);
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







