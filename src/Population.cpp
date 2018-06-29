/*
 * =====================================================================================
 *
 *       Filename:  Population.cpp
 *
 *    Description:  Class representing a population
 *
 *        Version:  1.0
 *        Created:  06/13/20178 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */


#include "Population.h"
#include "helpers.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <numeric>
#include <algorithm>
#include <assert.h>     /* assert */
#include <map>
#include <math.h>
#include <random>
#include <iterator>
#include <cmath>

Population::Population(){}

Population::Population(double r): rate(r) {}

Population::~Population(){}

double Population::getRate(){
	return rate;
}

void Population::setRate(double r){
	rate = r;
}

void Population::addUpdate(double prob, Update update){
	if(updates.size() != 0 && update.get().size() != updates[0].get().size()){
		throw("Intended update vector does not have same dimension as previous update vector");
	} else {
		probs.push_back(prob);
		updates.push_back(update);
		normalizeProbs();
	}
}

void Population::normalizeProbs(){
	norm_probs = normalize(probs);
}

std::vector<int> Population::getUpdate(double prob){
	if(prob < 0 || prob > 1){
		throw("Invalid probability given to population::getUpdate");
	}

	for(size_t i = 0; i < norm_probs.size(); i++){
		if(prob < norm_probs[i])
			return updates[i].get();
	}
}

void Population::printUpdates(){
	if(probs.size() != norm_probs.size() || norm_probs.size() != updates.size()){
		throw("Differing number of updates and update probilities");
	}

	for(size_t i = 0; i < updates.size(); i++){
		std::vector<int> up_i = updates[i].get();
		std::cout << "Update " << i << "; prob = " << norm_probs[i] << " {" << up_i[0];
		for(size_t j = 1; j < up_i.size(); j++){
			std::cout << ", " << up_i[j];
		}
		std::cout << "}" << std::endl;
	}
}
