#include "miscVDep.h"
using namespace std;

extern double xrylThresh;
extern double amphThresh;
extern PythonTools pyTools;
extern DataTools dataTools;
extern AtomTools atomTools;
extern SphHarmTools sphTools;	


bool atomType(Atom *a) {
	return a->type == "1" || a->type == "2"; 
}

bool type1(Atom *a) {
	return a->type == "1";
}

bool type2(Atom *a) {
	return a->type == "2";
}

bool alwaysTrue(Atom *a) {
	return true;
}

bool isXryl(Atom *a) {
	return a->qs[6] > xrylThresh;
}
bool isAmph(Atom *a) {
	return a->qs[6] < amphThresh;
}

bool isTrans(Atom *a) {
	return a->qs[6] >= amphThresh && a->qs[6] <= xrylThresh;
}

bool isAmphTrans(Atom *a) {
	return a->qs[6] <= xrylThresh;
}

bool diffType(Atom *a, Atom *n) {
	return n->type.compare(a->type) != 0;
}


void eraseExt(map<int, vector<string> > &fnsMap) {
	
	for (map<int, vector<string> >::iterator it=fnsMap.begin(); it!=fnsMap.end(); it++) {
		vector<string> &fns = it->second;
		cout << fns.size() << " files in group " << it->first << endl;
		for (unsigned int i=0; i<fns.size(); i++) {
			int pos = fns[i].find(".atom");
			fns[i].erase(pos, fns[i].size());
		}
	}
}

void calcExpStdev(experiment &e, function<double (run&)> fn, function<void (experiment&, double)> setVal) {
	double avg = 0;
	for (unsigned int i=0; i<e.runs.size(); i++) {
		avg += fn(e.runs[i]);
	}
	avg /= e.runs.size();
	double var = 0;
	for (unsigned int i=0; i<e.runs.size(); i++) {
		double dev = fn(e.runs[i]) - avg;
		var += dev * dev;
	}
	var /= e.runs.size();
	double stdev = sqrt(var);
	setVal(e, stdev);

}



void setAgVal(experiment& e, function<double (run&)> getVal, function<void (experiment&, double)> setExpVal, function<void (run&, double)> setRunVal) {
	double sumFrac;
	for (unsigned int i=0; i<e.runs.size(); i++) {
		run &r = e.runs[i];
		double frac = getVal(r);
		setRunVal(r, frac);
		sumFrac += frac;
	}
	setExpVal(e, sumFrac / e.runs.size());
}

void makeAtomGroup(experiment &e, run &r, bool (*testF) (Atom*), function<void (run&, vector<Atom *>&)> setGrp, function<void (experiment&, vector<Atom *>&)> insertAtoms, function<void (run&, double)> setFrac) {
	vector<Atom *> as = atomTools.test(&r.atomsInBounds, testF);
	setGrp(r, as);
	insertAtoms(e, as);
//	setFrac(r, as.size() / (double) r.atomsInBounds.size());
}

void assignQGroups(map<int, experiment> &exps) {
	for (map<int, experiment>::iterator it = exps.begin(); it!=exps.end(); it++) {
		experiment &e = it->second;
		for (unsigned int i=0; i<e.runs.size(); i++) {
			run &r = e.runs[i];
			MAKEATOMGROUP(e, r, isXryl, atomsXryl, fracXryl);
			MAKEATOMGROUP(e, r, isAmph, atomsAmph, fracAmph);
			MAKEATOMGROUP(e, r, isTrans, atomsTrans, fracTrans);
			MAKEATOMGROUP(e, r, isAmphTrans, atomsAmphTrans, fracAmphTrans);
		}
		SETAGVAL(e, fracXryl, r, r.atomsXryl.size() / (double) r.atomsInBounds.size()); 
		SETAGVAL(e, fracAmph, r, r.atomsAmph.size() / (double) r.atomsInBounds.size()); 
		SETAGVAL(e, fracTrans, r, r.atomsTrans.size() / (double) r.atomsInBounds.size()); 
		SETAGVAL(e, fracAmphTrans, r, r.atomsAmphTrans.size() / (double) r.atomsInBounds.size()); 
		/*e.fracXryl = e.atomsXryl.size() / total;
		e.fracAmph = e.atomsAmph.size() / total;
		e.fracTrans = e.atomsTrans.size() / total;
*/
		e.qXrylAvg = atomTools.avgQ(&e.atomsXryl, 6);
		e.qAmphAvg = atomTools.avgQ(&e.atomsAmph, 6);
		e.qTransAvg = atomTools.avgQ(&e.atomsTrans, 6);

		e.atomsAmphTrans.insert(e.atomsAmphTrans.end(), e.atomsAmph.begin(), e.atomsAmph.end());
		e.atomsAmphTrans.insert(e.atomsAmphTrans.end(), e.atomsTrans.begin(), e.atomsTrans.end());
		CALCEXPSTDEV(e, fracXryl, fracXrylStdev);
		CALCEXPSTDEV(e, fracAmph, fracAmphStdev);
		CALCEXPSTDEV(e, fracTrans, fracTransStdev);
		CALCEXPSTDEV(e, fracAmphTrans, fracAmphTransStdev);
	}

}

void populateQ6s(map<int, experiment> &exps) { //vector<vector<double> > &xss, vector<vector<double> > &yss, map<int, vector<string> > &fnsBase, map<int, vector<vector<Atom *> > > &atomsInBounds) {
	for (map<int, experiment>::iterator it=exps.begin(); it!=exps.end(); it++) {
		experiment &e = exps[it->first];
		vector<double> q6s; 
		for (unsigned int i=0; i<e.runs.size(); i++) {
			run &r = e.runs[i];
			for (unsigned int j=0; j<r.atomsInBounds.size(); j++) {
				q6s.push_back(r.atomsInBounds[j]->qs[6]);
			}
		}

		e.xssQ6Hist = vector<double>();
		e.yssQ6Hist = vector<double>();
		dataTools.hist(&q6s, &e.xssQ6Hist, &e.yssQ6Hist, .02, true);

	}
}

void makeTitles(map<int, experiment> &exps) {

	for (map<int, experiment>::iterator it=exps.begin(); it!=exps.end(); it++) {
		char buf[50];
		sprintf(buf, "%d", it->first);
		string title = "Tau = " + string(buf);
		it->second.title = title;
	}
}

void makeTitles(map<int, experiment> &exps, string prefix) {

	for (map<int, experiment>::iterator it=exps.begin(); it!=exps.end(); it++) {
		char buf[50];
		sprintf(buf, "%d", it->first);
		string title = prefix + string(buf);
		it->second.title = title;
	}
}


void populateMapMembers(map<int, experiment> &toPopulate, map<int, experiment> &src) {
	for (auto it=src.begin(); it!=src.end(); it++) {
		experiment e;
		toPopulate[it->first] = e;
	}
}

void processGrps(map<int, vector<string> > &fnsBase, map<int, experiment> &exps) {
	for (map<int, vector<string> >::iterator it=fnsBase.begin(); it!=fnsBase.end(); it++) {
		vector<string> *fnGroup = &it->second;
		experiment exp;
		exps[it->first] = exp;
		experiment &e = exps[it->first];
		cout << "starting new group " << it->first <<endl;
		e.id = it->first;	
		if (e.id == 1) {
			e.idStr = "V dep";
		} else {
			stringstream ss;
			ss << e.id;
			e.idStr = ss.str();
		}
		for (unsigned int i=0; i<fnGroup->size(); i++) {
			string fnBase = (*fnGroup)[i];
			run runny;
			e.runs.push_back(runny);
			run &r = e.runs[e.runs.size()-1];
			vector<Snap> snps = pyTools.loadDump(fnBase + ".atom", .003);	
			r.atoms = snps[0].atoms;
			r.grid = AtomGrid(&r.atoms, snps[0].bounds, 3, 3);
			r.grid.assignNeighbors(1.2, true, false);	
			for (unsigned int j=0; j<r.atoms.size(); j++) {
				r.atoms[j]->qs[6] = sphTools.calcQl(r.atoms[j], 6);
			}

			map<string, vector<double> > data = pyTools.grabData(fnBase + ".out", 1);
			r.eng = data["Eng"][0];
			
		}
		
	}
	cout << "loaded" << endl;
	for (map<int, experiment>::iterator it=exps.begin(); it!=exps.end(); it++) {
		experiment &e = it->second;
		for (unsigned int i=0; i<e.runs.size(); i++) {
			run &r = e.runs[i];
			vector<Atom *> isOfType = atomTools.test(&r.atoms, atomType);
			Bounds b = atomTools.findBounds(&isOfType);
			b.resize("y", .7, .5);
			r.atomsInBounds = atomTools.inBounds(&isOfType, b);
			e.type12.insert(e.type12.end(), isOfType.begin(), isOfType.end());	
			e.atomsInBounds.insert(e.atomsInBounds.end(), r.atomsInBounds.begin(), r.atomsInBounds.end());
		}
	}
}

void assignXrylGroups(map<int, experiment> &exps) {
	for (map<int, experiment>::iterator it=exps.begin(); it!=exps.end(); it++) {
		experiment &e = it->second;
//use floor > .45 . That would get most atoms
		for (unsigned int i=0; i<e.runs.size(); i++) {
			run &r = e.runs[i];
			r.crystalGroups = atomTools.assignCrystalGroups(r.atomsInBounds, 6, xrylThresh);
			e.crystalGroups.insert(e.crystalGroups.end(), r.crystalGroups.begin(), r.crystalGroups.end());
		}
	}


}

void linInterp(vector<double> &xs, vector<double> &ys, double *slope, double *intercept) {
	if (xs.size() != ys.size()) {
		cerr << "Doing regerssion with xs.size() != ys.size()" << endl;
	}

	double numPts = xs.size();
	double avgX = accumulate(xs.begin(), xs.end(), 0.0) / numPts;
	double avgY = accumulate(ys.begin(), ys.end(), 0.0) / numPts;
	double numer = 0;
	double denom = 0;
	for (int i=0; i<numPts; i++) {
		numer += (xs[i] - avgX) * (ys[i] - avgY);
		denom += (xs[i] - avgX) * (xs[i] - avgX);
	}
	if (denom == 0) {
		cerr << "Denom is zero in interp << endl" << endl;;
	}
	(*slope) = numer / denom;
	(*intercept) = avgY - (*slope) * avgX;
}

void interpLiqs(vector<map<int, experiment> > &srcs, map<int, experiment> interped) {
	
}


void convertToPys(map<int, experiment> &exps, PyObject *xssPy, PyObject *yssPy, PyObject *titlesPy) {
	for (map<int, experiment>::iterator it = exps.begin(); it!=exps.end(); it++) {
		experiment &e = it->second;
		PyObject *xsPy = pyTools.doubleVecToPy(&e.xssQ6Hist);
		PyObject *ysPy = pyTools.doubleVecToPy(&e.yssQ6Hist);
		PyObject *titlePy = PyString_FromString(e.title.c_str());
		PyList_Append(xssPy, xsPy);
		PyList_Append(yssPy, ysPy);
		PyList_Append(titlesPy, titlePy);

	}
}

void convertToPysXrylSize(map<int, experiment> &exps, PyObject *xss, PyObject *yss, PyObject *titlesPy) {
	for (map<int, experiment>::iterator it = exps.begin(); it!=exps.end(); it++) {
		experiment &e = it->second;

		PyObject *xsPy = pyTools.doubleVecToPy(&e.xsXrylSizeHist);
		PyObject *ysPy = pyTools.doubleVecToPy(&e.ysXrylSizeHist);
		PyObject *titlePy = PyString_FromString(e.title.c_str());
		PyList_Append(titlesPy, titlePy);
		PyList_Append(xss, xsPy);
		PyList_Append(yss, ysPy);

	}
}

void assignXrylHist(map<int, experiment> &exps) {

	for (map<int, experiment>::iterator it = exps.begin(); it!=exps.end(); it++) {
		experiment &e = it->second;
		vector<int> sizes;
		for (unsigned int i=0; i<e.crystalGroups.size(); i++) {
			sizes.push_back(e.crystalGroups[i]->size());
		}
		e.xsXrylSizeHist = vector<double>();
		e.ysXrylSizeHist = vector<double>();
		dataTools.intHist(&sizes, &e.xsXrylSizeHist, &e.ysXrylSizeHist, true);
	}

}

void push_backExps(map<int, experiment> &exps, vector<experiment *> &all) {

	for (map<int, experiment>::iterator it = exps.begin(); it!=exps.end(); it++) {
		all.push_back(&it->second);
	}
}




