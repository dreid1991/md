#ifndef ATOM_H
#define ATOM_H

#include "Vector.h"
#include <string>
#include <iostream>
#include <vector>
#include <complex>
#include "stamped_item.h"
#include "timestamped_pos_def.h"
#include <map>
using namespace std;
class Atom {
	public:
		Atom(double x_, double y_, double z_, int type_, int id_) : pos(Vector(x_, y_, z_)) {
			type = type_;
			id = id_;
			r = 0;
			crystalGroup = NULL;
			m = 1;
		}
		Atom(double x_, double y_, double z_, int type_, int id_, double m) : pos(Vector(x_, y_, z_)) {
			type = type_;
			id = id_;
			r = 0;
			crystalGroup = NULL;
			m = m;
		}
		vector<Atom *> ghosts;
		vector<Vector> ghostOffsets;
		vector<Atom*> neighbors;
		Vector accel;
		Vector accelLast;
		Vector vel;
		Vector pos;	
		Vector posWhenGrid;
		double r;
		double id;
		double m;
		vector<timestamped_pos> posHist;	
		string type;
		map<int, vector<complex<double> > > qms;
		map<int, double> qs;
		vector<Atom *> *crystalGroup;
		Vector vToNeighbor(int i) {
			return Vector(neighbors[i]->pos.x + neighborOffsets[i].x - pos.x, neighbors[i]->pos.y + neighborOffsets[i].y - pos.y);
		}
		void print() {
			cout << "\nPrintingAtom\n";
			cout << "id: ";
			cout << id;
			cout << "\ntype: ";
			cout << type;
			cout << "\nx: ";
			cout << pos.x;
			cout << "\ny: ";
			cout << pos.y;
			cout << "\n";
		}
};

#endif
