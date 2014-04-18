

void runStepInit(vector<Atom *> atoms, Box box, double ljCut, double neighCut, double epsilon, double sigma) {
	unsigned int ii = atoms.size();
	double ljCutSqr = ljCut * ljCut;
	for (unsigned int i=0; i<ii-1; i++) {
		for (unsigned int j=i+1; j<ii; j++) {
			Atom *a = atoms[i];
			Atom *b = atoms[j];
			Vector aPos = a->pos;
			Vector bPos = b->pos;
			if (fabs(aPos.x - bPos.x) > box.trace.x / 2) {
				if (aPos.x > bPos.x) {
					aPos.x -= box.trace.x;
				} else {
					aPos.x += box.trace.x;
				}
			}
			if (fabs(aPos.y - bPos.y) > box.trace.y / 2) {
				if (aPos.y > bPos.y) {
					aPos.y -= box.trace.y;
				} else {
					aPos.y += box.trace.y;
				}
			}
			if (fabs(aPos.z - bPos.z) > box.trace.z / 2) {
				if (aPos.z > bPos.z) {
					aPos.z -= box.trace.z;
				} else {
					aPos.z += box.trace.z;
				}
			}
			double distSqr = aPos.distSqr(bPos);
		
			if (distSqr < ljCutSqr) {
				dist = sqrt(distSqr);
				fMag = 4 * epsilon * (6 * pow(sigma, 6)	/ pow(distSqr, 3.5) - 12 * pow(sigma, 12) / pow(distSqr, 6.5));
				Vector force = Vector((bPos.x - aPos.x) * fMag / dist, (bPos.y - aPos.y) * fMag / dist, (bPos.y - aPos.y) * fmag / dist);
				a->accel += force / a->m;
				b->accel -= force / a->m;
			}
		}
	}
}
