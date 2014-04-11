#include "AtomTools.h"


vector<double> AtomTools::crystalGroupsTilt(vector<vector<Atom *> *> &crystalGroups) {
	vector<double> weights;
	for (unsigned int i=0; i<crystalGroups.size(); i++) {
		weights.push_back(crystalGroupTilt((*crystalGroups[i])));
	}
	return weights;
}

double AtomTools::crystalGroupTilt(vector<Atom *> &crystalGroup) { 
	vector<double> xs = LISTMAP(Atom *, double, a, crystalGroup, a->pos.x);
	vector<double> ys = LISTMAP(Atom *, double, a, crystalGroup, a->pos.y);
	double slope = dataTools.linSlope(xs, ys);
	double mag = sqrt(slope * slope + 1);
	Vector slopeUV = Vector(1 / mag, slope / mag);
	Vector perpSlopeUV = slopeUV;
	perpSlopeUV.rotatecw();
	double xWeight = 0;
	double yWeight = 0;
	for (unsigned int j=0; j<crystalGroup.size(); j++) {
		Atom &a = (*crystalGroup[j]);
		double xWeightAtomUnsqr = a.pos.dot(slopeUV) * a.m;
		double yWeightAtomUnsqr = a.pos.dot(perpSlopeUV) * a.m;
		xWeight += xWeightAtomUnsqr * xWeightAtomUnsqr;
		yWeight += yWeightAtomUnsqr * yWeightAtomUnsqr;
	}
	return xWeight / yWeight;	
}


vector<vector<int> > AtomTools::voronoiPolygonVertexIdxss(PyObject *voronoiData) {
	PyObject *polygonsPy = PyObject_GetAttrString(voronoiData, "regions");
	vector<vector<int> > polygonVertexIdxss = LISTMAPFROMPY(vector<int>, polygon, polygonsPy, 
			LISTMAPFROMPY(int, ptIdx, polygon, (int) PyInt_AsLong(ptIdx)) //PARSE THIS
			); 
//scipy's voronoi function returns voronoi object.  One of its attrs is the list of all involved vertices
//Another is the vertex indices that make of each polygon (list of list of ints).  Is -1 if goes off grid
	for (auto it = polygonVertexIdxss.end()-1; it>=polygonVertexIdxss.begin(); it--) {
		vector<int> &polygonVertexIdxs = (*it);
		if (polygonVertexIdxs.size() == 0) {
			polygonVertexIdxss.erase(it);
		}
		for (unsigned int j=0; j<polygonVertexIdxs.size(); j++) {
			if (polygonVertexIdxs[j] < 0) {
				polygonVertexIdxss.erase(it); //was boundary region, can't count points on polygon
				break;
			}
		}
	}
	/*
	for (unsigned int i=0; i<polygonVertexIdxss.size(); i++) {
		cout << "begin poly" << endl;
		for (unsigned int j=0; j<polygonVertexIdxss[i].size(); j++) {
			cout << polygonVertexIdxss[i][j] << endl;
		}
	}
	*/
	return polygonVertexIdxss;
}

vector<int> AtomTools::voronoiTypeCount(vector<Atom *> &atoms) {
	PyObject *voronoiData = pyTools.callVoronoi(atoms);
	vector<vector<int> > polyVertIdxss = voronoiPolygonVertexIdxss(voronoiData);
	return voronoiTypeCount(polyVertIdxss);
}

vector<int> AtomTools::voronoiTypeCount(vector<vector<int> > &polygonVertexIdxss) {

	vector<int> counts;
	auto expandCounts = [&] (const int x) {
		counts.reserve(x);
		for (int i=(signed int)counts.size(); i<=x; i++) {
			counts.push_back(0);
		}
	};
	for (unsigned int i=0; i<polygonVertexIdxss.size(); i++) {
		unsigned int polyType = polygonVertexIdxss[i].size();
		if (polyType > counts.size()) {
			expandCounts(polyType);
		}
		counts[polyType]++;
	}
	return counts;
}

vector<Polygon> AtomTools::voronoiPolygons(PyObject *voronoiData, vector<vector<int> > &vertexIdxss) {


	PyObject *verticesNumpy = PyObject_GetAttrString(voronoiData, "vertices"); 
	PyObject *verticesPy = pyTools.convertNPArray(verticesNumpy, true);
	vector<Vector> allVertices = LISTMAPFROMPY(Vector, vertex, verticesPy, 
				Vector(
					PyFloat_AsDouble(PyList_GetItem(vertex, 0)), 
					PyFloat_AsDouble(PyList_GetItem(vertex, 1))
				)
			);
	

	for (unsigned int i=0; i<allVertices.size(); i++) {
		allVertices[i].print();
	}
	vector<Polygon> polygons;
	for (unsigned int i=0; i<vertexIdxss.size(); i++) {
		polygons.push_back(Polygon());
		Polygon &poly = polygons[polygons.size() - 1];
		vector<int> &vertexIdxs = vertexIdxss[i];
		for (unsigned int j=0; j<vertexIdxs.size(); j++) {
			poly.push_back(allVertices[vertexIdxs[j]]);
		}
	}
	return polygons;
}




Bounds AtomTools::findBounds(vector<Atom *> *atoms) {
	double xlo = INFINITY;
	double xhi = -INFINITY;
	double ylo = INFINITY;
	double yhi = -INFINITY;
	for (unsigned int i=0; i<atoms->size(); i++) {
		xlo = fmin((*atoms)[i]->pos.x, xlo);
		xhi = fmax((*atoms)[i]->pos.x, xhi);
		ylo = fmin((*atoms)[i]->pos.y, ylo);
		yhi = fmax((*atoms)[i]->pos.y, yhi);
	}
	return Bounds(xlo, xhi, ylo, yhi);
}


void AtomTools::removeAtom(Atom *a, vector<Atom *> *as) {
	vector<Atom *>::iterator it = find(as->begin(), as->end(), a);
	if (it != as->end()) {
		as->erase(it);
	} 
}

void AtomTools::loopAtomHist(vector<timestamped_pos> &tsps, unsigned int curTspIdx, Bounds &b, bool loopX, bool loopY) {
	timestamped_pos &cur = tsps[curTspIdx];
	timestamped_pos &prev = tsps[curTspIdx-1]; 
	Vector offset;
	double xSpan = b.span("x");
	double ySpan = b.span("y");
	if (loopX) {
		if (cur.item.x > prev.item.x) {
			offset.x = xSpan;
		} else {
			offset.x = -xSpan;
		}
	}
	if (loopY) {
		if (cur.item.y > prev.item.y) {
			offset.y = ySpan;

		} else {
			offset.y = -ySpan;
		}
	}
	unsigned int startIdx = 0;
	unsigned int endIdx = curTspIdx;
	if (curTspIdx > ((double) tsps.size()) / 2) { //if over halfway through, move the later atoms instead of first
		startIdx = curTspIdx;
		endIdx = tsps.size();
		offset.x *= -1;
		offset.y *= -1;
	}
	for (unsigned int i=startIdx; i < endIdx; i++) {
		tsps[i].item += offset; 
	}
}

void AtomTools::loopPositionHists(vector<Atom *> &lastAtoms, Bounds b) {
	double xSpan = b.span("x");
	double ySpan = b.span("y");
	for (unsigned int i=0; i<lastAtoms.size(); i++) {
		vector<timestamped_pos> &tsps = lastAtoms[i]->posHist;
		for (unsigned int j=1; j<tsps.size(); j++) {
			timestamped_pos &cur = tsps[j];
			timestamped_pos &prev = tsps[j-1];
			bool loopX = fabs(cur.item.x - prev.item.x) > .5 * xSpan;
			bool loopY = fabs(cur.item.y - prev.item.y) > .5 * ySpan;
			if (loopX || loopY) {
				loopAtomHist(tsps, j, b, loopX, loopY); 
			}
		}
	}
}



//atoms are sorted by id.  See loadDump in PythonTools
void AtomTools::assignAtomPositionHist(vector<Snap> &snaps) {
	Snap &last = snaps[snaps.size()-1];
	for (unsigned int i=0; i<snaps.size(); i++) {
		Snap &s = snaps[i];
		vector<Atom *> &atoms = s.atoms;
		double timestamp = s.timestamp;
		for (unsigned int j=0; j<atoms.size(); j++) {
			timestamped_pos tsp;
			tsp.item = atoms[j]->pos;
			tsp.stamp = timestamp;
			last.atoms[j]->posHist.push_back(tsp);	
		}
	}
	Bounds assumedBounds = snaps[snaps.size()-1].bounds;
	loopPositionHists(last.atoms, assumedBounds);
}

void AtomTools::calcMobilityProfile(MobilityProfile &profile, vector<Atom *> &atoms, int stepIdx) {
	for (unsigned int i=0; i<atoms.size(); i++) {
		Atom *a = atoms[i];
		profile.addMovement(a->posHist[stepIdx], a->posHist[stepIdx+1]);

	}
	profile.calcProfile();
}

//will return vector of timestamped mobility profiles
vector<stamped_item<double, MobilityProfile> > AtomTools::calcMobilityProfiles(vector<Atom *> &atoms, double ylo, double yhi, double dy) {
	//all atoms' position histories taken from same snaps, so timestamps are all the same
	//that being said, this next move is kind of icky
	vector<stamped_item<double, MobilityProfile> > stamped_profiles;
	if (atoms.size() == 0) {
		return stamped_profiles;
	}
	int numProfiles = atoms[0]->posHist.size()-1;
	for (int i=0; i<numProfiles; i++) {
		double time = (atoms[0]->posHist[i].stamp + atoms[0]->posHist[i+1].stamp) / 2;
		stamped_item<double, MobilityProfile> stampedProfile;
		stampedProfile.stamp = time;
		stampedProfile.item = MobilityProfile(ylo, yhi, dy);
		calcMobilityProfile(stampedProfile.item, atoms, i);
		
		stamped_profiles.push_back(stampedProfile);
	}
	return stamped_profiles;		
}

vector<pair<double, displacement_datum> > AtomTools::calcDisplacementData(vector<Atom *> &atoms) {
	vector<pair<double, displacement_datum> > disps;
	auto getDispAtTime = [&] (const double time) {
		int idx = 0;
		for (auto it = disps.begin(); it !=disps.end(); it++) {
			if ((*it).first == time) {

				return &disps[it - disps.begin()];
			} 
			else if ((*it).first > time) {
				disps.insert(it, pair<double, displacement_datum> (time, displacement_datum_init()));
				return &disps[idx];
			}
			idx ++;
		}

		disps.push_back(pair<double, displacement_datum> (time, displacement_datum_init()));
		return &disps[disps.size()-1];
			
	};
	for (unsigned int i=0; i<atoms.size(); i++) {
		vector<timestamped_pos > &tsps = atoms[i]->posHist;
		for (unsigned int j=0, jj=tsps.size()-1; j<jj; j++) {
			for (unsigned int k=j+1, kk=jj+1; k<kk; k++) {
				timestamped_pos &earlier = tsps[j];
				timestamped_pos &later = tsps[k];

				double dt = later.stamp - earlier.stamp;	

				double rndDt = roundDouble(dt, 5);
				if (rndDt< 0) {
					cout << dt << " " << rndDt << endl;

				}
				pair<double, displacement_datum> *dispPair = getDispAtTime(rndDt);
				dispPair->second.displacements.push_back(later.item.distSqr(earlier.item));
			}
		}
	}
	avgStdevDisplacement(disps);
	return disps;
}


void AtomTools::avgStdevDisplacement(vector<pair<double, displacement_datum> > &disps) {
	for (unsigned int j=0; j<disps.size(); j++) {
		displacement_datum &dispDat = disps[j].second;
		vector<double> &dispVals = dispDat.displacements;
		double avg = 0;
		double stdev = 0;
		for (unsigned int i=0, ii=dispVals.size(); i<ii; i++) {
			avg += dispVals[i];
		}
		avg /= dispVals.size();
		for (unsigned int i=0, ii=dispVals.size(); i<ii; i++) {
			stdev += (dispVals[i] - avg) * (dispVals[i] - avg);
		}
		stdev = sqrt(stdev / dispVals.size());
		dispDat.stderror = stdev / dispVals.size();
		dispDat.avg = avg;
		dispDat.stdev = stdev;
	}
}

void AtomTools::growSeed(Atom *seed, vector<Atom *> *seeds, int qType, double qMin) {
	vector<Atom *> *myGrp = seed->crystalGroup;
	//a growing group should never run into an existing group, since any growing group will gobble up all the atoms it can reach with q>qMin
	//THEREFORE, I DECLARE that if an atom I run into has a crystal group assigned, it must be the same as the one I am growing
	for (unsigned int i=0; i<seed->neighbors.size(); i++) {
		Atom *neighbor = seed->neighbors[i];
		if (neighbor->qs[qType] >= qMin && neighbor->crystalGroup == NULL) {
			myGrp->push_back(neighbor);
			neighbor->crystalGroup = myGrp;
			removeAtom(neighbor, seeds);
			growSeed(neighbor, seeds, qType, qMin); 
		}	
	}
}

vector<vector<Atom *> *> AtomTools::assignCrystalGroups(vector<Atom *> seeds, int qType, double qMin) {
	vector<vector<Atom *> *> grps;
	for (int i= (int) seeds.size()-1; i>=0; i--) {
		if (seeds[i]->qs[qType] < qMin) {
			removeAtom(seeds[i], &seeds);
		}
	}

	while (seeds.size()) {
		Atom *seed = seeds[0];
		removeAtom(seed, &seeds);
		vector<Atom *> *newGrp = new vector<Atom *>;
		grps.push_back(newGrp);
		seed->crystalGroup = newGrp;
		newGrp->push_back(seed);
		growSeed(seed, &seeds, qType, qMin);
	}
	return grps;
}


vector<Atom *> AtomTools::inBounds(vector<Atom*> *atoms, Bounds b) {
	vector<Atom *> passed;
	for (unsigned int i=0; i<atoms->size(); i++) {
		if (b.atomInBounds((*atoms)[i])) {
			passed.push_back((*atoms)[i]);
		}
	}
	return passed;
}

double AtomTools::avgQ(vector<Atom *> *as, int q) {
	double sum = 0;
	for (unsigned int i=0; i<as->size(); i++) {
		sum += (*as)[i]->qs[q];
	}
	return sum / as->size();
}


vector<Atom *> AtomTools::selectByGrpSize(vector<vector<Atom *> *> &crystalGroups, int sizeMin, int sizeMax) {
	vector<Atom *> selected;
	if (sizeMax == -1) {
		sizeMax = INT_MAX;
	}
	for (unsigned int i=0; i<crystalGroups.size(); i++) {
		vector<Atom *> *grp = crystalGroups[i];
		int size = grp->size();
		if (size >= sizeMin && size <= sizeMax) {
			selected.insert(selected.end(), grp->begin(), grp->end());
		}
	}
	return selected;
}

double AtomTools::testNeighborsFrac(vector<Atom *> *atoms, bool (*testFSelf) (Atom *), bool (*testF) (Atom *, Atom *)) {
	int totalNeighs = 0;
	int totalPassed = 0;
	for (unsigned int i=0; i<atoms->size(); i++) {
		if (testFSelf((*atoms)[i])) {
			vector<Atom *> &neighs = (*atoms)[i]->neighbors;
			for (unsigned int j=0; j<neighs.size(); j++) {
				if (testF((*atoms)[i], neighs[j])) {
					totalPassed++;
				}
			}
			totalNeighs += neighs.size();
		}
	}
	return totalPassed / (double) totalNeighs;
}

int AtomTools::count(const vector<Atom *> *atoms, bool (*testF) (Atom *)) {
	int passed = 0;
	for (unsigned int i=0; i<atoms->size(); i++) {
		if (testF((*atoms)[i])) {
			passed++;
		}
	}
	return passed;
}


vector<Atom *> AtomTools::test(const vector<Atom *> *atoms, bool (*testF) (Atom *)) {
	vector<Atom *> passed;
	for (unsigned int i=0; i<atoms->size(); i++) {
		if (testF((*atoms)[i])) {
			passed.push_back((*atoms)[i]);
		}
	}
	return passed;
}
