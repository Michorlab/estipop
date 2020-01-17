#include "test_custom2.h"
// line above must point to header (.h) file with same filename as current file
double rate(double time, void* p){
	
	double ret = 0.2 * exp(-.6*time);
	
	// perform some maniupations...
	
	return ret;
}


