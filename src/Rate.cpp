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

Rate::Rate(){}

Rate::Rate(double (*f)(double, void*)){
	funct.function = f;
	rate_homog = maximizeFunc(funct, 0, 1000, 1000);
}

Rate::~Rate(){}

double Rate::integrateFunct(double a, double b){
	gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(1000);
	
	double result, error;
	size_t neval;

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
	  std::cout<<"Result " << result << " +/- " << error << " from " << neval << " evaluations" <<
		std::endl;
		return result;
	}
}

