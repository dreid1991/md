#include "displacement_datum.h"

displacement_datum displacement_datum_init() {
	displacement_datum d;
	d.avg = 0;
	d.stdev = 0;
	d.stderror = 0;
	return d;
}
