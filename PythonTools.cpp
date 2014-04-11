#include "PythonTools.h"

using namespace std;


string appendSlash(string s) {
	if (s[s.size()-1] != '/') {
		s += "/";
	}
	return s;
}

string cutLineEnd(string s) {
	while (s.size() > 0 && (s[s.size()-1] == '\n' || s[s.size()-1] == '\r')) {
		s.erase(s.size()-1, 1);
	}
	return s;
}

string cutEndSpaces(string s) {
	while (s.size() > 0 && s[s.size()-1] == ' ') {
		s.erase(s.size()-1, 1);
	}
	return s;
}

vector<double> getNums(string s) {
	vector<double> nums;
	istringstream iss(s);
	
	double num;
	while (iss >> num) {
		nums.push_back(num);
	} 
	return nums;
}

void appendData(vector<vector<double> > *cols, vector<double> nums) {
	if (nums.size() != cols->size()) {
		cout << "Num of data columns changed" << endl;
	}
	for (unsigned int i=0; i<nums.size(); i++) {
		(*cols)[i].push_back(nums[i]);
	}
}

void PythonTools::setDump(PyObject *dumpFunc_) {
	dumpFunc = dumpFunc_;
}

void PythonTools::init() {
//TERRIBLE: IMPLIED THAT WE ADDED PATH FOR EXTRAS
	PyObject *extras = PyImport_ImportModule("extras");
	plt = PyImport_ImportModule("matplotlib.pyplot");
	re = PyImport_ImportModule("re");
	reCompile = PyObject_GetAttrString(re, "compile");
	reSearch = PyObject_GetAttrString(re, "search");
	plot = PyObject_GetAttrString(plt, "plot");
	show = PyObject_GetAttrString(plt, "show");
	savefig = PyObject_GetAttrString(plt, "savefig");
	scatter = PyObject_GetAttrString(plt, "scatter");
	hist = PyObject_GetAttrString(plt, "hist");

	isNone = PyObject_GetAttrString(extras, "isNone");

	scipy = PyImport_ImportModule("scipy");
	spatial = PyImport_ImportModule("scipy.spatial");
	voronoi = PyObject_GetAttrString(spatial, "Voronoi");
}

void PythonTools::setPlt(PyObject *plt_) {
	//lol deprecated
}

PyObject *PythonTools::callVoronoi(vector<Atom *> &atoms) {
	auto makePtList = [] (const Atom *a) {
		PyObject *myList = PyList_New(2);
		PyList_SetItem(myList, 0, PyFloat_FromDouble(a->pos.x));
		PyList_SetItem(myList, 1, PyFloat_FromDouble(a->pos.y));
		return myList;
	};


	PyObject *atomPts = LISTMAPPY(Atom *, a, atoms, makePtList(a));
	return PyObject_CallFunction(voronoi, (char *)"O", atomPts);
}
void PythonTools::setGrab(PyObject *grabData_) {
	grabDataMod = grabData_;
	Py_INCREF(grabDataMod);
}

vector<map<string, double> > PythonTools::grabOutComments(string folderName, string fileName, vector<string> commentsToRead) {
	folderName = appendSlash(folderName);	
	FILE *fr = fopen((folderName + fileName).c_str(), "rt");
	char line[150];

	vector<vector<double> > vals;
	
	for (unsigned int i=0; i<commentsToRead.size(); i++) {
		vals.push_back(vector<double>());
	}
	while (fgets(line, 150, fr) != NULL) {
		string lineStr = string(line);
		for (unsigned int i=0; i<commentsToRead.size(); i++) {
			string comment = commentsToRead[i];
			if (lineStr.find(comment) != string::npos) {
				lineStr.erase(lineStr.begin(), lineStr.begin() + comment.size());
				double value = atof(lineStr.c_str());
				vals[i].push_back(value);
			}
		}
	}
	vector<map<string, double> > comments;
	auto buildComment = [&] (const int &idx) {
		map<string, double> comment;
		for (unsigned int j=0; j<vals.size(); j++) {
			comment[commentsToRead[j]] = vals[j][idx];
		}
		return comment;
	};
	unsigned int numEntries = vals[0].size();
	for (unsigned int i=0; i<vals.size()-1; i++) {
		if (vals[i].size() != vals[i+1].size()) {
			cout << "Reading comments from file " << fileName << " in folder " << folderName << " and comment lengths do not match" << endl;
			cout << commentsToRead[i] << " is of length " << vals[i].size() << " and " << commentsToRead[i+1] << " is of length " << vals[i+1].size() << endl;
		}
	}
	for (unsigned int i=0; i<numEntries; i++) {
		comments.push_back(buildComment(i));
	}
	fclose(fr);
	return comments;
}


map<int, PyObject *> PythonTools::compileRes(map<int, string> &reGroups) {
	map<int, PyObject *> pyRes;
	for (auto it=reGroups.begin(); it!=reGroups.end(); it++) {
		int grpId = it->first;
		string reExpr = it->second;
		PyObject *pyRe = PyObject_CallFunction(reCompile, (char *) "s", reExpr.c_str());
		pyRes[grpId] = pyRe;
	}
	return pyRes;

}

map<int, vector<string> > PythonTools::grabFileNames(string folderName, map<int, string> reGroups) {
	cout << "folder name is " << folderName << endl;
	folderName = appendSlash(folderName);
	cout.flush();
	map<int, vector<string > > data;
	for (auto it=reGroups.begin(); it!=reGroups.end(); it++) {
		data[it->first] = vector<string>();
	}

	DIR *dp;
	struct dirent *dirp;
	if ((dp = opendir(folderName.c_str())) == NULL) {
		cerr << "can't open folder " << folderName << endl;
	}
	while ((dirp = readdir(dp)) != NULL) {
		string fn = string(dirp->d_name);

		for (map<int, vector<string> >::iterator it=data.begin(); it!= data.end(); it++) {
			PyObject *result = PyObject_CallFunction(reSearch, (char *) "ss", reGroups[it->first].c_str(), fn.c_str());
			PyObject *noneRes = PyObject_CallFunction(isNone, (char *) "O", result);
			//segfaults and confusion await all who touch this code
			//doing re.comile just... doesn't work
			int isNone = PyLong_AsLong(noneRes);
			if (!isNone) {
				data[it->first].push_back(folderName + fn);
			}

		}
	}


	return data;
}

bool atomIdLess(Atom *a, Atom *b) {
	return a->id < b->id;
}

vector<Snap> PythonTools::loadDump(string fn, double timestep) {
	vector<Snap> snaps;
	FILE *fr = fopen(fn.c_str(), "rt");
	cout << "opened file " << endl;
	char line[150];
	vector<string> lines;
	Snap *s;
	double xlo, xhi, ylo, yhi;
	
	while(fgets(line, 150, fr) != NULL) {
		Snap ss;
		snaps.push_back(ss);
		s = &snaps[snaps.size()-1];
		fgets(line, 150, fr);
		s->timestamp = atof(line) * timestep;
		fgets(line, 150, fr);
		fgets(line, 150, fr);
		int numAtoms = atoi(line);
		fgets(line, 150, fr);
		fgets(line, 150, fr);
		vector<double> xs = getNums(cutLineEnd(string(line)));
		xlo = xs[0];
		xhi = xs[1];
		fgets(line, 150, fr);
		vector<double> ys = getNums(cutLineEnd(string(line)));
		ylo = ys[0];
		yhi = ys[1];	
		fgets(line, 150, fr);
		fgets(line, 150, fr);
		for (int i=0; i<numAtoms; i++) {
			fgets(line, 150, fr);
			vector<double> atomBits = getNums(cutLineEnd(string(line)));
			double x = xlo + (xhi - xlo) * atomBits[2];
			double y = ylo + (yhi - ylo) * atomBits[3];
			ostringstream strs;
			strs << atomBits[1];
			
			string type = strs.str();
			int id = (int) atomBits[0];
			Atom *a = new Atom(x, y, type, id);
			s->atoms.push_back(a);
			
		}
		sort(s->atoms.begin(), s->atoms.end(), atomIdLess);
		s->box= Box(xlo, xhi, ylo, yhi);
	}
	fclose(fr);
	return snaps;


}

map<int, string> PythonTools::readDescData(string folder, string entry) {
	folder = appendSlash(folder);	
	FILE *fr = fopen((folder + "runDesc.txt").c_str(), "rt");
	char line[150];
	int runNum = 0;
	map<int, string> vals;
	while (entry[entry.size()-1] == ' ') {
		entry.erase(entry.end()-1, entry.end());
	}
	while (fgets(line, 150, fr) != NULL) {
		string lineStr = string(line);
		if (lineStr.find("RUN ") != string::npos) {
			lineStr.erase(lineStr.begin(), lineStr.begin() + 4);
			runNum = atoi(lineStr.c_str());
		} else if (lineStr.find(entry) != string::npos) {
			int eraseTo = lineStr.find(':') + 2;
			lineStr.erase(lineStr.begin(), lineStr.begin() + eraseTo);
			vals[runNum] = cutEndSpaces(cutLineEnd(lineStr));
		}
	}
	fclose(fr);
	return vals;

}
vector<string> PythonTools::grabVDepByDesc(vector<string> folders, string descItem, string val) {
	vector<string> validFns;
	for (unsigned int i=0; i<folders.size(); i++) {
		string folder = appendSlash(folders[i]);
		map<int, string> descData = readDescData(folder, descItem);
		vector<int> validRuns;
		for (map<int, string>::iterator it = descData.begin(); it != descData.end(); it++) {
			if (val.compare(it->second) == 0) {
				validRuns.push_back(it->first);
			}
		}
		map<int, string> searchMap;
		for (unsigned int i=0; i<validRuns.size(); i++) {
			string testNum = intToStr(validRuns[i]);
			string reExpr = "minend" + testNum + "-[0-9]+\\.atom"; 	
			searchMap[i] = reExpr;
		}

		map<int, vector<string> > fns = grabFileNames(folder + "postProc", searchMap);
		for (map<int, vector<string> >::iterator it = fns.begin(); it!=fns.end(); it++) {
			validFns.insert(validFns.end(), it->second.begin(), it->second.end());
		}
	}
	return validFns;

}


void PythonTools::scatterAtoms(vector<Atom *> &as, PyObject *dict) {
	vector<double> xs;
	vector<double> ys;
	xs.reserve(as.size());
	ys.reserve(as.size());
	for (unsigned int i=0; i<as.size(); i++) {
		xs.push_back(as[i]->pos.x);
		ys.push_back(as[i]->pos.y);
	}
	PyObject *xsPy = doubleVecToPy(&xs);
	PyObject *ysPy = doubleVecToPy(&ys);
	PyObject *args = Py_BuildValue("OO", xsPy, ysPy);
	PyObject_Call(scatter, args, dict);
	Py_DECREF(xsPy);
	Py_DECREF(ysPy);



}

bool PythonTools::allNums(string str) {	
	if (str.size() == 0) {
		return false;
	}
	for (unsigned int i=0; i<str.size(); i++) {
		if (!((str[i]>='0' && str[i] <= '9') || str[i] == 'e' || str[i] == 'E' || str[i] == '-' || str[i] == '.' || str[i] == ' ')) {
			return false;
		}
	}
	return true;
}

int PythonTools::countDataCols(vector<string> *lines) {
	int numCols = 0;
	for (int i = (int) lines->size()-1; i>=0; i--) {
		if (allNums((*lines)[i])) {
			istringstream iss((*lines)[i]);
			do {
				string sub;
				iss >> sub;
				if (!iss) {
					break;
				}
				numCols++;
			} while (iss);
			return numCols;
		}
	}
	cerr << "CAN'T COUNT COLUMNS LOL";
	return -1;
}
map<string, vector<double> > PythonTools::grabData(string fn, long numGrps) {

	FILE *fr = fopen(fn.c_str(), "rt");
	char line[150];
	vector<string> lines;
	fgets(line, 150, fr);
	while(fgets(line, 150, fr) != NULL) {
		lines.push_back(cutLineEnd(string(line)));
	}
	int numDataCols = countDataCols(&lines);
	vector<vector<double> > cols;
	for (int i=0; i<numDataCols; i++) {
		cols.push_back(vector<double>());
	}
	bool inData = false;
	bool inMinimization = false;
	int curNumGrps = 0;
	string titleLine = "";
	for (int i=(int) lines.size()-1; i>=0; i--) {
		bool isAllNums = allNums(lines[i]);
		if (isAllNums && inData && !inMinimization) {
			appendData(&cols, getNums(lines[i]));
			curNumGrps++;
		} else if (!isAllNums && inData) { //proxy for in header
			titleLine = lines[i];
			if (numGrps == curNumGrps) {
				break;
			}
			inData = false;
		} else if (lines[i].find("Loop time of") != string::npos && !inMinimization) {
			inData = true;
		} else if (lines[i].find("Minimization stats:") != string::npos) {
			inMinimization = true;
		} else if (lines[i].find("Setting up minimization ...") != string::npos) {
			inMinimization = false;
		}
	}
	for (unsigned int i=0; i<cols.size(); i++) {
		reverse(cols[i].begin(), cols[i].end());
	}
	vector<string> titles;
	istringstream iss(titleLine);
	string title;
	while (iss >> title) {
		titles.push_back(title);
	}
	map<string, vector<double> > data;
	for (unsigned int i=0; i<titles.size(); i++) {
		data[titles[i]] = cols[i];
	}
	fclose(fr);
	return data;
}


	
string PythonTools::vectorToStr(vector<double> *v) {
	string s = "[";
	for (unsigned int i=0; i<v->size()-1; i++) {
		ostringstream strs;
		strs<<(*v)[i];
		string str = strs.str();
		s += str + ",";
	}
	ostringstream strs;
	strs<<(*v)[v->size()-1];
	string str = strs.str();
	s += str +  "]";
	return s;
}

PyObject *PythonTools::doubleVecToPy(vector<double> *xs) {
	PyObject *vec = PyList_New(xs->size());
	for (unsigned int i=0; i<xs->size(); i++) {
		PyObject *x = PyFloat_FromDouble((*xs)[i]);
		PyList_SetItem(vec, i, x);
//		Py_DECREF(x); //if you decref x, the value in the list goes to NULL or something.  So don't do that, I guess
	}
	return vec;
}



PyObject *PythonTools::mapIntDoubleToPy(map<int, double> &m) {
	PyObject *dict = PyDict_New();
	for (map<int, double>::iterator it = m.begin(); it!=m.end(); it++) {
		PyObject *key = PyInt_FromLong(it->first);
		PyObject *val = PyFloat_FromDouble(it->second);
		PyDict_SetItem(dict, key, val);
	}
	return dict;
}

vector<string> PythonTools::PyToVecString(PyObject *p) {
	int size = PyList_Size(p);
	vector<string> c;
	for (int i=0; i<size; i++) {
		PyObject *a = PyList_GetItem(p, i);
		c.push_back(PyString_AsString(a));
		Py_DECREF(a);
	}
	return c;

}

vector<double> PythonTools::PyToVecDouble(PyObject *p) {
	int size = PyList_Size(p);
	vector<double> c;
	for (int i=0; i<size; i++) {
		PyObject *a = PyList_GetItem(p, i);
		c.push_back(PyFloat_AsDouble(a));
		Py_DECREF(a);
	}
	return c;
}

PyObject * PythonTools::convertNPArray(PyObject *numpyArray, bool decrefNumpy) {
	PyObject *tolist = PyObject_GetAttrString(numpyArray, "tolist");
	PyObject *listPy = PyObject_CallFunction(tolist, (char *)""); 
	Py_DECREF(tolist);
	if (decrefNumpy) {
		Py_DECREF(numpyArray);
	}
	return listPy;
}
