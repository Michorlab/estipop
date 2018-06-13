/*
 * =====================================================================================
 *
 *       Filename:  State.cpp
 *
 *    Description:  Class representing system state
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


#include "State.h"
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
#include <gsl/gsl_randist.h>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;

State::State(){
	
}

State::State(std::vector<int> s){
	vec = s;
	
}

State::~State(){
	
}

void State::print(){
	for(size_t i = 0; i < vec.size(); i++){
		std::cout << vec[i] << "\t";
	}
	std::cout << std::endl;
}

void State::updateState(std::vector<int> update){
	if(update.size() != vec.size()){
		throw std::invalid_argument("States do not have same dimension");
	}
	
	for(size_t i = 0; i < update.size(); i++){
		vec[i] += update[i];
	}
}

void State::addPopulation(Population p){
	pops.push_back(p);
}

double State::getNextTime(){
	double overall = 0.0;
	for(size_t i = 0; i < pops.size(); i++){
		overall += vec[i] * pops[i].getRate();
	}
	std::cout << overall << std::endl;
	return(gsl_ran_exponential(rng, 1 / overall));
}

int State::choosePop(){
	std::vector<double> rates;
	for(size_t i  = 0; i < pops.size(); i++){
		rates.push_back(vec[i] * pops[i].getRate());
	}
	
	int choice = choose(rates);
	return choice;
}