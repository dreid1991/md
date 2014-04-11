#ifndef GRID_H
#define GRID_H

#include <vector>
using namespace std;
template <class T>
class Grid {
	public:
		vector<vector<T> > cols;
		Grid() {
			
		}
		int nx;
		int ny;
		Grid(int nx_, int ny_) {
			nx = nx_;
			ny = ny_;
			for (int x=0; x<nx; x++) {
				vector<T> col;
				for (int y=0; y<ny; y++) {
					T sqr;
					col.push_back(sqr);
				}
				cols.push_back(col);
		
			}
		}
		unsigned int size() {
			return cols.size();
		}
		vector<T> &operator[](int x) {
			return cols[x];
		}
};

#endif
