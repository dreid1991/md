
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
Vector AtomGrid::sqrPosition(vector <Atom*> *sqr) {
	for (unsigned int i=0; i<grid.size(); i++) {
		for (unsigned int j=0; j<grid[i].size(); j++) {
			for (unsigned int k=0; k<grid[i][j].size(); k++) {
				if (&grid[i][j][k] == sqr) {
					return Vector(box.xlo + dx * i, box.ylo + dy * j, box.zlo + dz * k);
				}
			}
		}
	}
	return Vector(-1, -1, -1);
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







