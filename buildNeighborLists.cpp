#include "buildNeighborLists.h"


void buildNeighborLists(vector<Atom *> &atoms, Box box, int nProcX, int nProcY, int nProcZ, double ljCut, double neighCut) {
	double dxProc = box.trace.x / nProcX;
	double dyProc = box.trace.y / nProcY;
	double dzProc = box.trace.z / nProcZ;
	AtomGrid cpuGrid = AtomGrid(atoms, box, dxProc, dyProc, dzProc, false);
	int numProc = nProcX * nProcY * nProcZ;
	vector<AtomGrid> hiResGrids;
	hiResGrids.reserve(numProc);
	for (unsigned int x=0; x<cpuGrid.nx; x++) {
		for (unsigned int y=0; y<cpuGrid.ny; y++) {
			for (unsigned int z=0; z<cpuGrid.nz; z++) {
				Box boxBounds = cpuGrid.getSquareBox(x, y, z);
				hiResGrids.push_back(AtomGrid(boxBounds, 2 * neighCut, 2 * neighCut, 2 * neighCut));
				hiResGrids[hiResGrids.size()-1].populateAssignNeigh(cpuGrid[x][y][z], neighCut);
			}
		}
	}
	//now implement ghosts/looping
}




/*


void assignNeighFromSquares(Grid<vector<Atom *> > &grid, Atom *a, int xIdx, int yIdx, int zIdx) {
	for (int dx=-1; dx<=1; dx++) {
		for (int dy=-1; dy<=1; dy++) {
			for (int dz=-1; dz<=1; dz++) {
				int x = xIdx + dx - (int) grid.nx * floor((xIdx + dx) / grid.nx);
				int y = yIdx + dy - (int) grid.ny * floor((yIdx + dy) / grid.ny);
				int z = zIdx + dz - (int) grid.nz * floor((zIdx + dz) / grid.nz);
				vector<Atom *> &gridSqr = grid[x][y][z];
				a->neighbors.erase(a->neighbors.begin(), a->neighbors.end());
				a->neighbors.insert(a->neighbors.end(), gridSqr.begin(), gridSqr.end());
			}
		}
	}
}

void AtomTools::assignNeighbors(vector<Atom *> atoms, double neighCut, Box box) {
	int nx = ceil(box.trace.x / (2 * neighCut));
	int ny = ceil(box.trace.y / (2 * neighCut));
	int nz = ceil(box.trace.z / (2 * neighCut));
	double gridDim = 2 * neighCut;
	Grid<vector<Atom *> > grid = Grid<vector<Atom *> >(nx, ny, nz);
	for (unsigned int i=0; i<atoms.size(); i++) {
		Atom *a = atoms[i];
		box.loopIntoBox(a);
		int xIdx = (a->pos.x - box.xlo) / gridDim;
		int yIdx = (a->pos.y - box.ylo) / gridDim;
		int zIdx = (a->pos.z - box.zlo) / gridDim;
		assignNeighFromSquares(grid, a, xIdx, yIdx, zIdx);	
		grid[xIdx][yIdx][zIdx].push_back(a);
		a->posWhenGrid = a->pos;
	}
}
*/
