#include "buildNeighborLists.h"



void assignFaceNeighbors(vector<vector<Atom *> &> assignToSqrs, vector<vector<Atom *> &> assignFromSqrs, int nMajor, int nMinor, Vector offset, double cutSqr) {
	for (unsigned int i=0; i<assignToSqrs.size(); i++) {
		vector<Atom *> &sqr = assignToSqrs[i];
		int major = i / nMajor;
		int minor = i - nMajor * major;
		for (int fromMajor = fmax(major-1, 0), int fromMajorMax = fmin(major+2, nMajor); fromMajor < fromMajorMax; fromMajor++) {
			for (int fromMinor = fmax(minor-1, 0), int fromMinorMax = fmin(minor+2, nMinor); fromMinor < fromMinorMax; fromMinor++) {
				vector<Atom *> &toSqr = assignToSqrs[nMajor * major + minor];
				vector<Atom *> &fromSqr = assignFromSqrs[fromMajor * nMajor + fromMinor];

				for (unsigned int a=0, aa = toSqr.size(); a < aa; a++) {
					for (unsigned int b=0; bb = fromSqr.size(); b<bb; b++) {
						Atom *atomA = toSqr[a];
						Atom *atomB = fromSqr[b];
						if (atomA->pos.distSqr(atomB->pos + offset) < cutSqr) {
							atomA->ghosts.push_back(atomB);
							atomA->ghostOffsets.push_back(offset);

						}
					}
					
				}
			}
			
		}
	}
}

void assignGhostsInDirection(AtomGrid &cpuGrid, vector<AtomGrid> &hiresGrids, int x, int y, int z, int dx, int dy, int dz, int checkingX, int checkingY, int checkingZ, Vector offset, double cutSqr) {
	AtomGrid &assignTo = hiresGrids[x * cpuGrid.ny * cpuGrid.nz + y * cpuGrid.nz + z]; //because high res's were flattened
	AtomGrid &assignFrom = hiresGrids[checkingX * cpuGrid.ny * cpuGrid.nz + checkingY * cpuGrid.nz + checkingZ];
	if (fabs(dx) + fabs(dy) + fabs(dz) == 1) { //faces touching
		tuple<vector<vector<Atom *> &>, int, int> assignToSqrs = assignTo.selectFace(assignTo.faceId(dx, dy, dz));
		tuple<vector<vector<Atom *> &>, int, int> assignFromSqrs = assignFrom.selectFace(assignTo.faceId(-dx, -dy, -dz));
		//you were RIGHT HERE.  need to implement selectFace and get it to return nMajor, nMinor.  Those are the grid dimensions of the face and are used to map the 1D vector onto 2D grid idxs
		assignFaceNeighbors(assignToSqrs, assignFromSqrs, offset, cutSqr);
	} else if (fabs(dx) + fabs(dy) + fabs(dz) == 2) { //edges touching

	} else { //corners touching

	}
}

void assignGhosts(AtomGrid &cpuGrid, vector<AtomGrid> &hiresGrids, int x, int y, int z, double cut) {
	//need to check for ghosts in all neighboring cpu grids
	for (int dx=-1; dx<=1; dx++) {
		for (int dy=-1; dy<=1; dy++) {
			for (int dz=-1; dz<=1; dz++) {
				if (!(dx==0 && dy==0 && dz==0)) {
					int xOffset = -cpuGrid.nx * floor((x + dx) / cpuGrid.nx);
					int yOffset = -cpuGrid.ny * floor((y + dy) / cpuGrid.ny);
					int zOffset = -cpuGrid.nz * floor((z + dz) / cpuGrid.nz);
					Vector offset = Vector(cpuGrid.box.dx * xOffset, cpuGrid.box.dy + yOffset, cpuGrid.box.dz);
					
					int checkingX = x + xOffset;
					int checkingY = x + yOffset;
					int checkingZ = x + zOffset;
					assignGhostsInDirection(cpuGrid, hiresGrids, x, y, z, dx, dy, dz, checkingX, checkingY, checkingZ, offset, cut * cut);
				}
			}
		}
	}
}

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
				hiResGrids.push_back(AtomGrid(boxBounds, 2 * neighCut, 2 * neighCut, 2 * neighCut, true));
				hiResGrids[hiResGrids.size()-1].populateAssignNeigh(cpuGrid[x][y][z], neighCut);
			}
		}
	}
	for (unsigned int x=0; x<cpuGrid.nx; x++) {
		for (unsigned int y=0; y<cpuGrid.ny; y++) {
			for (unsigned int z=0; z<cpuGrid.nz; z++) {
				assignGhosts(cpuGrid, hiResGrids, x, y, z, neighCut);
			}
		}
	}
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
