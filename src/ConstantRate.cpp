/*
 * =====================================================================================
 *
 *       Filename:  ConstantRate.cpp
 *
 *    Description:  Time-dependent rate functions
 *
 *        Version:  1.0
 *        Created:  09/13/2018 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */


#include "ConstantRate.h"
#include "helpers.h"

#include <iostream>
#include <fstream>
#include <gsl/gsl_randist.h>

extern gsl_rng* rng;

double constantRate(double x, void* p){
	constant_params &params= *reinterpret_cast<constant_params *>(p);
	return (params.rate);
}

double linearRate(double x, void* p){
	linear_params &params= *reinterpret_cast<linear_params *>(p);
	return (params.intercept + x * params.slope);
}


ConstantRate::ConstantRate(double r){
	params.rate = r;
	
	funct.function = &constantRate;
	funct.params = reinterpret_cast<void *>(&params);
	
	rate_homog = r;
}

ConstantRate::~ConstantRate() {}

double ConstantRate::operator()(double time){
	return params.rate;
}

double ConstantRate::operator()(double begin, double end){
	return (end - begin) * params.rate;
}

LinearRate::LinearRate(double i, double s){
	params.intercept = i;
	params.slope = s;

	funct.function = &linearRate;
	funct.params = reinterpret_cast<void *>(&params);

	rate_homog = maximizeFunc(funct, 0, 1000, 1000);
}

LinearRate::~LinearRate() {}

double LinearRate::operator()(double time){
	return params.intercept + params.slope * time;
}

double LinearRate::operator()(double begin, double end){
	double result, error;
	size_t neval;

	const double xlow=begin;
	const double xhigh=end;
	const double epsabs=1e-4;
	const double epsrel=1e-4;

	int code=gsl_integration_qng(&funct,
								xlow,
								xhigh,
								epsabs,
								epsrel,
								&result,
								&error,
								&neval);
	 if(code){
		std::cerr<<"There was a problem with integration: code " << code
             <<std::endl;
	} else {
		tot_error += error;
		return result;
	}

	return -1;
}

