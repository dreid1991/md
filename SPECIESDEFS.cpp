#include <map>
#include "SpeciesDef.h"
#include "PARMS.h"
using namespace std;
void assertSpeciesDefs() {
	PARMS::SPECIESDEFS.insert(make_pair(1, SpeciesDef(1, 2, 1, 1)));

//	SPECIESDEFS[1] = a
}
