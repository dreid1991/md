#ifndef GRID_H
#define GRID_H

#include <vector>
using namespace std;
template <class T>
class Grid {
	public:
		vector<vector<vector<T> > > xslices;
		Grid() {
			
		}
		int nx;
		int ny;
		int nz;
		Grid(int nx_, int ny_, int nz_) {
			xslices.reserve(nx);
			nx = nx_;
			ny = ny_;
			nz = nz_;
			for (int x=0; x<nx; x++) {
				vector<vector<T> > xslice;
				xslice.reserve(ny);
				for (int y=0; y<ny; y++) {
					vector<T> ycol;
					ycol.reserve(nz);
					for (int z=0; z<nz; z++) {
						T sqr;
						ycol.push_back(sqr);
					}
					xslices.push_back(ycol);
				}
				xslices.push_back(xslice);
		
			}
		}
		unsigned int size() {
			return xslices.size();
		}
		vector<vector<T> > &operator[](int x) {
			return xslices[x];
		}
};

#endif
