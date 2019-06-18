/*
 * =====================================================================================
 *
 *       Filename:  Update.cpp
 *
 *    Description:  Class representing a system update
 *
 *        Version:  1.0
 *        Created:  06/28/2018 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */


#include "Update.h"
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

Update::Update(){}

Update::Update(std::vector<int> f) : fixed(f){
	is_random = false;
}

Update::Update(bool rand, std::vector<double> o, std::string dist, std::vector<double> p) : is_random(rand), offspring_vec(o), offspring_dist(dist), params(p){

}

Update::~Update(){}

/*
std::vector<int> Update::get(){
	if(!is_random){
		return fixed;
	}

	for(size_t i = 0; i < random_indices.size(); i++){
		fixed[random_indices[i]] = gsl_ran_poisson(rng, 5) + 1;
	}

	return fixed;
}
*/

std::vector<int> Update::get(){
	if(!is_random){
		return fixed;
	}

	int offspring;

	if(offspring_dist == "poisson"){
		offspring = gsl_ran_poisson(rng, params[0]);
	} else {
		offspring = 2;
	}

	int K = offspring_vec.size();
	std::vector<unsigned int> ret(K);

	gsl_ran_multinomial(rng, K, offspring, offspring_vec.data(), ret.data());

	std::vector<int> ret_int(ret.begin(), ret.end());

	return ret_int;
}
