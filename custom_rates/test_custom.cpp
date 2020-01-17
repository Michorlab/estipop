#include "test_custom.h"
// line above must point to header (.h) file with same filename as current file
double rate(double time, void* p){
	
	double ret = 0.3 * exp(-.6*time);
	
	// perform some maniupations...
	
	return ret;
}


