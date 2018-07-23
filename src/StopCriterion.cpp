/*
 * =====================================================================================
 *
 *       Filename:  StopCrition.cpp
 *
 *    Description:  Class representing a stop to the system
 *
 *        Version:  1.0
 *        Created:  07/23/2018 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */

#include "StopCriterion.h"
#include "helpers.h"

#include <iostream>
#include <fstream>
//#include <iomanip>
//#include <stdlib.h>     /* srand, rand */
//#include <time.h>       /* time */
//#include <numeric>
//#include <algorithm>
//#include <assert.h>     /* assert */
//#include <map>
//#include <math.h>
//#include <random>
//#include <iterator>
//#include <cmath>
#include <gsl/gsl_randist.h>

extern gsl_rng* rng;

//

StopCriterion::StopCriterion(){}

StopCriterion::StopCriterion(std::vector<int> ind, std::string ineq, double v) : indices(ind), inequality(ineq), value(v){}

StopCriterion::~StopCriterion(){}

bool StopCriterion::check(std::vector<int> state){
	int checkVal = 0;
	for(size_t i = 0; i < indices.size(); i++){
		checkVal += state[indices[i]];
	}

	/*
	switch(inequality){
		case ">": {
			if(checkVal > value)
				return true;
			return false;
		} case ">=": {
			if(checkVal >= value)
				return true;
			return false;
		} case "<": {
			if(checkVal < value)
				return true;
			return false;
		} case "<=": {
			if(checkVal <= value)
				return true;
			return false;
		} case "=": {
			if(checkVal == value)
				return true;
			return false;
		}
	}*/

	if(inequality == ">"){
		if(checkVal > value)
			return true;
		return false;
	} else if (inequality == ">="){
		if(checkVal >= value)
			return true;
		return false;
	}else if (inequality == "<"){
		if(checkVal < value)
			return true;
		return false;
	}else if (inequality == "<="){
		if(checkVal <= value)
			return true;
		return false;
	}else if (inequality == "="){
		if(checkVal == value)
			return true;
		return false;
	}


	return false;
}


