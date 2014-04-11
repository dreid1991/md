#ifndef MOBILITY_PROFILE_H
#define MOBILITY_PROFILE_H
#include <vector>
#include "timestamped_pos_def.h"
#include "math.h"
#include <utility>
#include "Vector.h"
#include "stamped_item.h"
using namespace std;
struct mobility_slice {
	double y;
	double dy;
	double mobility;
};

class MobilityProfile {
	vector<pair<timestamped_pos, timestamped_pos> > movementBuffer;
	public:
		MobilityProfile(double ylo_, double yhi_, double dy_);
		MobilityProfile();
		vector<mobility_slice> profile;
		void addMovement(timestamped_pos &, timestamped_pos &);
		void calcProfile();
		double mobilityAt(double y);
		double ylo, yhi, dy;
		vector<pair <timestamped_pos, timestamped_pos> > movements; 
};


#endif

