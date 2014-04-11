#include "MobilityProfile.h"

MobilityProfile::MobilityProfile(double ylo_, double yhi_, double dy_) {
	ylo = ylo_;
	yhi = yhi_;
	dy = dy_;
}

MobilityProfile::MobilityProfile() {

}

void MobilityProfile::addMovement(timestamped_pos &a, timestamped_pos &b) {
	movements.push_back(pair<timestamped_pos, timestamped_pos> (a, b));

}

void MobilityProfile::calcProfile() {
	vector<double> weights;
	unsigned int numSlices = ceil((yhi - ylo) / (double) dy);
	for (unsigned int i=0; i<numSlices; i++) {
		mobility_slice m;
		m.y = ylo + i * dy;
		m.dy = dy;
		m.mobility = 0;
		profile.push_back(m);
		weights.push_back(0.);
	}

	for (unsigned int i=0; i<movements.size(); i++) {
		double dt =	fabs(movements[i].first.stamp - movements[i].second.stamp);
		Vector a = movements[i].first.item;
		Vector b = movements[i].second.item;
		//cout << "ya is " << a.y << endl;
		//cout << "yb is " << b.y << endl;
		double yloPair = fmin(a.y, b.y);
		double yhiPair = fmax(a.y, b.y);
		double mobility = a.distSqr(b) / dt;
		int idxlo = (yloPair - ylo) / dy;
		int idxhi = (yhiPair - ylo) / dy;
		double weight = 1 / (double) (idxhi - idxlo + 1);
		//cout << "calced weight is " << weight << endl;
		//cout << "idx lo is " << idxlo << " idx hi is " << idxhi << endl;
		//cout << "dt is " << dt << endl;
		for (int j=idxlo; j<=idxhi; j++) {
			weights[j] += weight;
			profile[j].mobility += mobility * weight;
			//cout << "weight of " << j << " is " << weights[j] << endl;
			//cout << "mobility is " << profile[j].mobility << endl;
		}
	}
	
	for (unsigned int i=0; i<numSlices; i++) {
		if (weights[i] > 0) {
			profile[i].mobility /= weights[i];
		}
	}

}
