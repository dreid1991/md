#ifndef FILESPEC_H
#define FILESPEC_H
#include <string>
using namespace std;
class FileSpec {
	public:
		string prefix;
		int testMin;
		int testMax;
		int runMin;
		int runMax;
		string suffix;
		FileSpec(string prefix_, int testMin_, int testMax_, int runMin_, int runMax_, string suffix_) {
			prefix = prefix_;
			testMin = testMin_;
			testMax = testMax_;
			runMin = runMin_;
			runMax = runMax_;
			suffix = suffix_;
		}
		FileSpec() {

		}
};

#endif
