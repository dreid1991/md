void updatePositions(vector<Atom *> &atoms, double dt) {
	for (unsigned int i=0, unsigned int ii = atoms.size(); i<ii; i++) {
		Atom *a = atoms[i];
		a->pos += dt * a->vel + (.5 * dt * dt) * a->accel;
	}
}


void updateVelocities(vector<Atom *> &atoms, double dt) {
	for (unsigned int i=0, unsigned int ii = atoms.size(); i<ii; i++) {
		Atom *a = atoms[i];
		a->vel += (.5 * dt) * (a->accel + a->accelLast);
	}
}

void runStepInit(vector<Atom *> atoms, Box box, double ljCut, double neighCut, double epsilon, double sigma, double dt) {
	updatePositions(atoms, dt);
	unsigned int ii = atoms.size();
	double ljCutSqr = ljCut * ljCut;
	for (unsigned int i=0; i<ii-1; i++) {
		Atom *a = atoms[i];
		for (unsigned int j=i+1; j<ii; j++) {
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
		a->accelLast = a->accel;
	}
	updateVelocities(atoms, dt);
}


void runSteps(vector<Atom *> atoms, Box box, double ljCut, double neighCut, double epsilon, double sigma, double dt) {
	updatePositions(atoms, dt);
	unsigned int ii = atoms.size();
	double ljCutSqr = ljCut * ljCut;
	for (unsigned int i=0; i<ii-1; i++) {
		a->accelLast = a->accel;
		Atom *a = atoms[i];
		for (unsigned int j=i+1; j<ii; j++) {
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
	updateVelocities(atoms, dt);
}
