#ifndef SPECIESDEFS_H
#define SPECIESDEFS_H
#include <iostream>
#include <map>
#include "SpeciesDef.h"
#include "PARMS.h"
using namespace std;
namespace asserts {
	void assertSpeciesDefs() {
		PARMS::SPECIESDEFS.insert(make_pair(1, SpeciesDef(1, 2, 1, 1)));
		SpeciesDef foo = PARMS::SPECIESDEFS.find(1)->second;
		cout << foo.m << endl;
	
	}
}

#endif
