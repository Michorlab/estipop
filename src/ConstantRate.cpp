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
extern gsl_integration_workspace *workspace;

double constantRate(double x, void* p){
	constant_params &params= *reinterpret_cast<constant_params *>(p);
	return (params.rate);
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

// LinearRate

double linearRate(double x, void* p){
	linear_params &params= *reinterpret_cast<linear_params *>(p);
	return (params.intercept + x * params.slope);
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
	return std::max(0.0, params.intercept + params.slope * time);
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

// SwitchRate

double switchRate(double x, void* p){
	switch_params &params= *reinterpret_cast<switch_params *>(p);
	if(x < params.tswitch)
		return params.pre;
	else
		return params.post;
}

SwitchRate::SwitchRate(double pre, double post, double t){
	params.pre = pre;
	params.post = post;
	params.tswitch = t;

	funct.function = &switchRate;
	funct.params = reinterpret_cast<void *>(&params);

	rate_homog = maximizeFunc(funct, 0, 1000, 1000);
}

SwitchRate::~SwitchRate() {}

double SwitchRate::operator()(double time){
	if(time < params.tswitch)
		return params.pre;
	else
		return params.post;
}

double SwitchRate::operator()(double begin, double end){
	double result, error;

	const double xlow=begin;
	const double xhigh=end;
	const double epsabs=1e-4;
	const double epsrel=1e-4;

	int code=gsl_integration_qags (&funct,
                                xlow,
                                xhigh,
                                epsabs,
                                epsrel,
                                1000,
                                workspace,
                                &result,
                                &error);
	 if(code){
		std::cerr<<"There was a problem with integration: code " << code
             <<std::endl;
	} else {
		tot_error += error;
		return result;
	}

	return -1;
}

// Pulse

double pulseRate(double x, void* p){
	pulse_params &params= *reinterpret_cast<pulse_params *>(p);

	double a = floor(x);
	double diff = x - a;

	double t = fmod((int)a, params.totPeriod);

	if(t + diff < params.lowPeriod){
		return params.low;
	} else {
		return params.high;
	}
}

PulseRate::PulseRate(double totPeriod, double lowPeriod, double low, double high){
	params.totPeriod = totPeriod;
	params.lowPeriod = lowPeriod;
	params.low = low;
	params.high = high;

	funct.function = &pulseRate;
	funct.params = reinterpret_cast<void *>(&params);

	rate_homog = maximizeFunc(funct, 0, 1000, 1000);
}

PulseRate::~PulseRate() {}

double PulseRate::operator()(double time){
	double a = floor(time);
	double diff = time - a;

	double t = fmod((int)a, params.totPeriod);

	if(t + diff < params.lowPeriod){
		return params.low;
	} else {
		return params.high;
	}
}

double PulseRate::operator()(double begin, double end){
	double result, error;
	//size_t neval;

	const double xlow=begin;
	const double xhigh=end;
	const double epsabs=1e-4;
	const double epsrel=1e-4;

	int code=gsl_integration_qags (&funct,
                                xlow,
                                xhigh,
                                epsabs,
                                epsrel,
                                1000,
                                workspace,
                                &result,
                                &error);
	 if(code){
		std::cerr<<"There was a problem with integration: code " << code
             <<std::endl;
	} else {
		tot_error += error;
		return result;
	}

	return -1;
}
