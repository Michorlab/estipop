/*
 * =====================================================================================
 *
 *       Filename:  ConstantRate.h
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

#pragma once
#include <string>
#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>

#include <gsl/gsl_integration.h>
//#include <map>

#include "Rate.h"

struct constant_params{
	double rate;
};

struct linear_params{
	double slope;
	double intercept;
};

class ConstantRate : public Rate {
public:

constant_params params;

ConstantRate(double r);

~ConstantRate();

virtual double operator()(double time);

virtual double operator()(double begin, double end);
	
};

class LinearRate : public Rate {
public:

linear_params params;

LinearRate(double i, double s);

~LinearRate();

virtual double operator()(double time);

virtual double operator()(double begin, double end);
	
};

