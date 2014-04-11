#include "DataTools.h"
#include <iostream>

string intToStr(int x) {
	stringstream ss;
	ss << x;
	return ss.str();
}

double DataTools::linSlope(vector<double> &xs, vector<double> &ys) {
	if (xs.size() != ys.size()) {
		cout << "Trying to do linear regression" << endl;
	}
	int n = xs.size();
	double avgX = accumulate(xs.begin(), xs.end(), 0.0) / n;
	double avgY = accumulate(ys.begin(), ys.end(), 0.0) / n;

	double numerator = 0.0;
	double denominator = 0.0;

	for (int i=0; i<n; i++) {
		numerator +=  (xs[i] - avgX) * (ys[i] - avgY);
		denominator += (xs[i] - avgX) * (xs[i] - avgX);
	}

	return numerator / denominator;
}
void DataTools::binAvg(map<int, map<string, vector<double> > > *data, string xPath, string yPath, vector<double> *xsDest, vector<double> *ysDest, int numBins, bool doMovingAvg, int span) {
	//xsDest, ysDest must be empty
	vector<double> mergedXs, mergedYs;
	vector<vector<double> > binXs, binYs;
	for (map<int, map<string, vector<double> > >::iterator i=data->begin(); i!= data->end(); i++) {
		vector<double> *xs = &i->second[xPath];
		vector<double> *ys = &i->second[yPath];
		appendList(&mergedXs, xs);
		appendList(&mergedYs, ys);
	}
	for (int i=0; i<numBins; i++) {
		vector<double> leVector;
		binXs.push_back(leVector);
		binYs.push_back(leVector);
	}
	double xMin = INFINITY;
	double xMax = -INFINITY;
	for (unsigned int i=0; i<mergedXs.size(); i++) {
		xMin = fmin(xMin, mergedXs[i]);
		xMax = fmax(xMax, mergedXs[i]);
	}
	double dxBin = (xMax - xMin) / (numBins - 1);
	for (unsigned int i=0; i<mergedXs.size(); i++) {
		int binIdx = floor((mergedXs[i] - xMin) / dxBin);
		binXs[binIdx].push_back(mergedXs[i]);
		binYs[binIdx].push_back(mergedYs[i]);
	}
	for (unsigned int i=0; i<binXs.size(); i++) {
		double xSum;
		double ySum;
		for (unsigned int j=0; j<binXs[i].size(); j++) {
			xSum += binXs[i][j];
			ySum += binYs[i][j];
		}
		if (binXs[i].size()) {
			(*xsDest).push_back(xSum / binXs[i].size());
			(*ysDest).push_back(ySum / binXs[i].size());
		}
	}
	if (doMovingAvg) {
		vector<double> avgdYs = movingAvg(xsDest, ysDest, span);
		for (unsigned int i=0; i<avgdYs.size(); i++) {
			(*ysDest)[i] = avgdYs[i];
		}
	}

}

void DataTools::appendList(vector<double> *xs, vector<double> *as) {
	xs->reserve(xs->size() + as->size());
	for (unsigned int i=0; i<as->size(); i++) {
		xs->push_back((*as)[i]);
	}
}

double DataTools::triangleAvg(vector<double> *xs, vector<double> *ys, int i, int span) {
	double sumVals, totalWeight;
	int width = fabs((*xs)[i+span] - (*xs)[i-span]);
	double centerWeight = 1 / ((double) span);
	double weightSlope = -centerWeight / width;
	for (int j=i-span; j<=i+span; j++) {
		double weight = fmax(0, centerWeight + weightSlope * fabs((*xs)[j] - (*xs)[i]));
		sumVals += weight * (*ys)[j];
		totalWeight += weight;
	}
	return sumVals / totalWeight;
}

vector<double> DataTools::movingAvg(vector<double> *xs, vector<double> *ys, int span) {
	if (!xs->size() || xs->size() != ys->size()) {
		cerr << "Tried to do moving avg on mismatched or zero data set";
	}
	vector<double> avgd;
	avgd.reserve(xs->size());
	for (int i=0; i<span; i++) {
		double sumVals, totalWeight;
		double width = fabs((*xs)[span + i] - (*xs)[0]);
		double centerWeight = 1 / ((double) span + 1);//this should probably be + i, not 1 
		double weightSlope = -centerWeight / width;
		for (int j=0; j<=span+i; j++) {
			double weight = fmax(0, centerWeight + weightSlope * fabs((*xs)[j] - (*xs)[i]));
			sumVals += weight * (*ys)[j];
			totalWeight += weight;		
		}
		avgd.push_back(sumVals/totalWeight);
	}
	for (unsigned int i=span; i<xs->size() - span; i++) {
		avgd.push_back(triangleAvg(xs, ys, i, span));
	}
	for (unsigned int i=xs->size()-span; i<xs->size(); i++) {
		double sumVals, totalWeight;
		double width = fabs((*xs)[xs->size() - 1] - (*xs)[i - span]);
		double centerWeight = 1 / ((double) span);
		double weightSlope = -centerWeight / width;
		for (unsigned int j=i-span; j<xs->size(); j++) {
			double weight = fmax(0, centerWeight + weightSlope * fabs((*xs)[j] - (*xs)[i]));
			sumVals += weight * (*ys)[j];
			totalWeight += weight;		
		}
		avgd.push_back(sumVals/totalWeight);
	}
	return avgd;
}

void DataTools::hist(vector<double> *events, vector<double> *xs, vector<double> *ys, double xSpan, bool normalize) {
	
	sort(events->begin(), events->end());
	double weightSlope = 1 / xSpan; //with center value = 1
	int eventFlr = 0;
	double xMin = (*events)[0];
	double xMax = (*events)[events->size()-1];

	int numPts = 300;
	double xStep = (xMax - xMin) / (numPts - 1);
	for (double x=xMin; x<=xMax; x+= xStep) {
		double yVal = 0;
		for (int eventIdx = eventFlr; (*events)[eventIdx] < x + xSpan && eventIdx < (int) events->size(); eventIdx++) {
			if ((*events)[eventIdx] >= x - xSpan && (*events)[eventIdx] <= x + xSpan) {
				yVal += fmax(0, 1 - weightSlope * fabs(x - (*events)[eventIdx]));
			}
		}
		ys->push_back(yVal);
		xs->push_back(x);
		while ((*events)[eventFlr+1] < x+xStep-xSpan) {
			eventFlr++;
		}
	}
	if (normalize) {
		double sumInt;
		for (unsigned int i=0; i<xs->size()-1; i++) {
			sumInt += .5 * ((*ys)[i] + (*ys)[i+1]) * ((*xs)[i+1] - (*xs)[i]); //trap integration
		}
		for (unsigned int i=0; i<ys->size(); i++) {
			(*ys)[i] /= sumInt;
		}
	}
	
	
}

void DataTools::histNoTri(vector<double> *events, vector<double> *xs, vector<double> *ys, double xSpan, bool normalize) {
	
	sort(events->begin(), events->end());
	int eventFlr = 0;
	double xMin = (*events)[0];
	double xMax = (*events)[events->size()-1];

	int numPts = 300;
	double xStep = (xMax - xMin) / (numPts - 1);
	for (double x=xMin; x<=xMax; x+= xStep) {
		double yVal = 0;
		for (int eventIdx = eventFlr; (*events)[eventIdx] < x + xSpan && eventIdx < (int) events->size(); eventIdx++) {
			if ((*events)[eventIdx] >= x - xSpan && (*events)[eventIdx] <= x + xSpan) {
				yVal += 1;
			}
		}
		ys->push_back(yVal);
		xs->push_back(x);
		while ((*events)[eventFlr+1] < x+xStep-xSpan) {
			eventFlr++;
		}
	}
	if (normalize) {
		double sumInt;
		for (unsigned int i=0; i<xs->size()-1; i++) {
			sumInt += .5 * ((*ys)[i] + (*ys)[i+1]) * ((*xs)[i+1] - (*xs)[i]); //trap integration
		}
		for (unsigned int i=0; i<ys->size(); i++) {
			(*ys)[i] /= sumInt;
		}
	}
	
	
}



void DataTools::intHist(vector<int> *events, vector<double> *xs, vector<double> *ys, bool normalize) {
	int xMin = 1000000;
	int xMax = -1000000;
	for (unsigned int i=0; i<events->size(); i++) {
		xMin = fmin(xMin, (*events)[i]);
		xMax = fmax(xMax, (*events)[i]);
	}

	cout << "max " << xMax << endl;
	int totalHeight = 0;
	for (int i=xMin; i<=xMax; i++) {
		xs->push_back(i);
		int height = 0;
		for (unsigned int j=0; j<events->size(); j++) {
			if ((*events)[j] == i) {
				height++;
			}
		}
		totalHeight += height;
		ys->push_back((double) height);
	}
	if (normalize) {
		for (unsigned int i=0; i<ys->size(); i++) {
			(*ys)[i] /= totalHeight;
		}
	}

}
