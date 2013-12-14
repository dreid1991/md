#include <iostream>
#include "dog.h"
using namespace std;


//class Dog {
//	public:
//		int numLegs;
//		Dog(int);
//};
//Dog::Dog(int numLegs_) {
//	numLegs = numLegs_;
//}

int main() {
	Dog pooch = Dog(5);
	cout << pooch.numLegs << endl;
	cout << "Hello World!" << endl;
	return 0;

}
