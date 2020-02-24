/*
 * =====================================================================================
 *
 *       Filename:  Rate.cpp
 *
 *    Description:  Time-dependent rate functions
 *
 *        Version:  1.0
 *        Created:  09/10/2018 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */


#include "Rate.h"
#include "helpers.h"

#include <iostream>
#include <fstream>
#include <gsl/gsl_randist.h>

extern gsl_rng* rng;
extern bool silent;

Rate::Rate(){};

Rate::Rate(double (*f)(double, void*)) : Rate(){
	funct.function = f;
	funct.params = 0; //reinterpret_cast<void *>(&params);
	rate_homog = maximizeFunc(funct, 0, 1, 1000);
}

Rate::~Rate(){}

double Rate::eval(double time){
	return std::max(0.0, GSL_FN_EVAL(&funct, time));
}

double Rate::operator()(double time){
	return std::max(0.0, GSL_FN_EVAL(&funct, time));
}

