#include "Vector.h"
#include <stdlib.h>
#include "Box.h"

using namespace std;
Box::Box(Vector &origin_, Vector &trace_): origin(origin_), trace(trace_) {
	this->origin = origin;
	this->trace = trace;
}

bool Box::inside(Vector &v) {
	return v.x > origin.x && v.x < origin.x + trace.x && v.y > origin.y && v.y < origin.y + trace.y && v.z > origin.z && v.z < origin.z + trace.z;
}

Vector Box::randV() {
	return Vector(origin.x + ((double) rand() / RAND_MAX) * trace.x, origin.y + ((double) rand() / RAND_MAX) * trace.y, origin.z + ((double) rand() / RAND_MAX) * trace.z);
}


