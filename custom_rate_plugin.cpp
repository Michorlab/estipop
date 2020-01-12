#include "custom_rate_plugin.h"
// line above must point to header (.h) file with same filename as current file
// [[Rcpp::export]]
double rate(double time, void* p){
	
	double ret = 0.1 * time;
	
	// perform some maniupations...
	
	return ret;
}


