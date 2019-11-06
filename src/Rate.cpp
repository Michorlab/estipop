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
extern gsl_integration_workspace* workspace;
extern bool silent;

Rate::Rate(){tot_error = 0;}

Rate::Rate(double (*f)(double, void*)) : Rate(){
	funct.function = f;
	funct.params = reinterpret_cast<void *>(&params);
	rate_homog = maximizeFunc(funct, 0, 1000, 1000);
}

Rate::~Rate(){}

double Rate::integrateFunct(double a, double b){
	//gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

	double result, error;
	//size_t neval;

	const double xlow=a;
	const double xhigh=b;
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
    gsl_integration_workspace_free(workspace);

	if(code)
	{
	  std::cerr<<"There was a problem with integration: code " << code
			   <<std::endl;
	}
	else
	{
	  //std::cout<<"Result " << result << " +/- " << error << " from " << neval << " evaluations" <<
	//	std::endl;
		tot_error += error;
	}
	return result;
}

double Rate::eval(double time){
	return std::max(0.0, GSL_FN_EVAL(&funct, time));
}

double Rate::operator()(double time){
	return std::max(0.0, GSL_FN_EVAL(&funct, time));
}

double Rate::operator()(double a, double b){
	//gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);

	double result, error;
	//size_t neval;

	const double xlow=a;
	const double xhigh=b;
	const double epsabs=1e-10;
	const double epsrel=1e-10;

	int code=gsl_integration_qags (&funct,
                                xlow,
                                xhigh,
                                epsabs,
                                epsrel,
                                1000,
                                workspace,
                                &result,
                                &error);
    //gsl_integration_workspace_free(workspace);

	if(code)
	{
	  std::cerr<<"There was a problem with integration: code " << code
			   <<std::endl;
	}
	else
	{
		//if(!silent)
	  //std::cout<<"Result " << result << " +/- " << error << " from " << neval << " evaluations" <<
		//std::endl;
		tot_error += error;
	}
	return result;

}
