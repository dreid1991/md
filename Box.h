#include "Vector.h"

#ifndef BOX_H
#define BOX_H

class Box{
	public:
		Box(Vector &, Vector &);
		Vector origin;
		Vector trace;
		bool inside(Vector &);
		Vector randV();
};

#endif
