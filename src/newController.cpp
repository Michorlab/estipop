/*
* =====================================================================================
*
*       Filename:  newController.cpp
*
*    Description:  Contains simulating and input file parsing functions
*
*        Version:  1.0
*        Created:  06/13/2018 14:57:27
*       Revision:  none
*       Compiler:  g++
*
*         Author:  Jeremy Ferlic (s), jferlic@g.harvard.edu
*   Organization:  Harvard University
*
* =====================================================================================
*/

// For plugin system
#ifndef USE_PRECOMPILED_HEADERS
#ifdef _WIN32
#include <direct.h>
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#endif

// Classes
#include "System.h"
#include "Update.h"
#include "StopCriterion.h"
#include "Rate.h"
#include "ConstantRate.h"

// Includes
#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <sstream>
#include <iterator>
#include <chrono>
#include <cstring>
#include <cmath>
#include <limits>
#include <gsl/gsl_randist.h>
#include <cmath>

// Rcpp headers
#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

// GSL random number generators
gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
double seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

bool silent = false;



// Helper Methods to read in data from text files
template
< typename T
, template<typename ELEM, typename ALLOC=std::allocator<ELEM> > class Container
>
std::ostream& operator<< (std::ostream& o, const Container<T>& container)
{
	typename Container<T>::const_iterator beg = container.begin();

	o << "["; // 1

	while(beg != container.end())
	{
		o << " " << *beg++; // 2
	}

	o << " ]"; // 3

	return o;
}

// trim from left
static inline std::string &ltrim(std::string &s)
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(),
			std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s)
{
	s.erase(std::find_if(s.rbegin(), s.rend(),
			std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s)
{
	return ltrim(rtrim(s));
}


// Read in a vector matrix
std::vector<std::vector<double> > fileToVectorMatrix(std::string name)
{
	std::vector<std::vector<double> > result;
	std::ifstream input (name);
	std::string lineData;

	while(getline(input, lineData))
	{
		double d;
		std::vector<double> row;
		std::stringstream lineStream(lineData);

		while (lineStream >> d)
		row.push_back(d);

		result.push_back(row);
	}

	return result;
}

template <typename T>
std::vector<T> flatten(const std::vector<std::vector<T> >& v)
{
	std::size_t total_size = 0;
	for (const auto& sub : v)
	total_size += sub.size(); // I wish there was a transform_accumulate
	std::vector<T> result;
	result.reserve(total_size);
	for (const auto& sub : v)
	result.insert(result.end(), sub.begin(), sub.end());
	return result;
}

// Read in a vector
std::vector<double> fileToVector(std::string name)
{
	std::vector<std::vector<double> > result;
	std::ifstream input (name);
	std::string lineData;

	while(getline(input, lineData))
	{
		double d;
		std::vector<double> row;
		std::stringstream lineStream(lineData);

		while (lineStream >> d)
		row.push_back(d);

		result.push_back(row);
	}

	return flatten(result);
}

// Read in a text file, output a vector of strings, one string per line, trimmed
std::vector<std::string> inputStringVector(std::string in)
{
	std::string s;
	std::vector<std::string> a;
	std::ifstream infile(in);
	if(infile.is_open())
	{
		while(getline(infile, s))
		{
			// allow comments
			std::string::iterator end_pos = std::remove(s.begin(), s.end(), ' ');
			s.erase(end_pos, s.end());
			a.push_back(trim(s));
		}
	}

	return(a);
}

std::vector<std::vector<std::string>> funct(std::string filename){
	std::ifstream csv(filename);
	std::string line;
	std::vector <std::vector<std::string>> items;

	if (csv.is_open()) {
		for (std::string row_line; std::getline(csv, row_line);)
		{
			items.emplace_back();
			std::istringstream row_stream(row_line);
			for(std::string column; std::getline(row_stream, column, ',');)
			items.back().push_back(column);
		}
	}
	else {
		std::cout << "Unable to open file";
	}

	return items;
}

//' gmbp3
//'
//' gmbp3
//'
//' @export
// [[Rcpp::export]]
double gmbp3(Rcpp::NumericVector observations, std::string file, Rcpp::NumericVector initial, Rcpp::List transitions, Rcpp::List stops, bool silence, SEXP seed = R_NilValue){
	double seedcpp;
	if(Rf_isNull(seed)){
		seedcpp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	} else {
		seedcpp = Rf_asReal(seed);
	}
    gsl_rng_set(rng, seedcpp);
	silent = silence;

	if(!silent) std::cout << "Starting process... " << std::endl;
	int nTrans = transitions.length();
	int nStops = stops.length();

	if(!silent) std::cout << "Initialization system..." << std::endl;
	// Initial population sizes
	std::vector<long int> init(initial.begin(), initial.end());

	// Initialize system
	System sys(init);

	/*
	// Create populations with lifetimes
	std::cout << "Making populations..." << std::endl;
	for(size_t i = 0; i < init.size(); i++){
	std::cout << "init[" << i << "]: " << init[i] << std::endl;
	Population* p = new Population(lifetimes[i]);
	sys.addPopulation(p);
	}
	*/

	// Add transitions
	if(!silent) std::cout << "Adding transitions..." << std::endl;

	// Iterate over transitions list
	for(int i = 0; i < nTrans; i++){

		// Get the ith transition
		Rcpp::List list_i = Rcpp::as<Rcpp::List>(transitions[i]);

		// Get popuation, is_random, probability
		int population = list_i[0];
		bool is_random = list_i[1];
		double rate = list_i[2];

		// If there is a random component, get which indices of the vector will be generated randomly
		if(is_random){
			// Get the offspring distribution of the Transition
			Rcpp::NumericVector oVec = Rcpp::as<Rcpp::NumericVector>(list_i[3]);
			std::vector<double> offspringVec (oVec.begin(), oVec.end());

			// Get the distribution
			std::string dist = list_i[4];

			// Get the offspring distribution of the Transition
			Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(list_i[5]);
			std::vector<double> offspringParams (p.begin(), p.end());

			// Add this transition
			sys.addUpdate(rate, population, Update(is_random, offspringVec, dist, offspringParams));
		} else{
			// Get the fixed portion of the Transition
			Rcpp::NumericVector fix = Rcpp::as<Rcpp::NumericVector>(list_i[3]);
			std::vector<int> fixed (fix.begin(), fix.end());

			// add update
			sys.addUpdate(rate, population, Update(fixed));
		}
	}

	// Add transitions
	if(!silent) std::cout << "Adding stopping criteria..." << std::endl;

	// Iterate over transitions list
	for(int i = 0; i < nStops; i++){

		// Get the ith transition
		Rcpp::List list_i = Rcpp::as<Rcpp::List>(stops[i]);

		// Get popuation, is_random, probability

		// Get the offspring distribution of the Transition
		Rcpp::NumericVector ind = Rcpp::as<Rcpp::NumericVector>(list_i[0]);
		std::vector<int> indices (ind.begin(), ind.end());

		std::string ineq = list_i[1];
		double value = list_i[2];

		sys.addStop(StopCriterion(indices, ineq, value));
	}

	//sys.print();

	// Observation times
	std::vector<double> obsTimes(observations.begin(), observations.end());

	// Simulate
	if(!silent) std::cout << "Simulating..." << std::endl;
	sys.simulate(obsTimes, file);
	if(!silent) std::cout << "Ending process..." << std::endl;

	return 0.0;
}


//' timeDepBranch
//'
//' timeDepBranch
//'
//' @export
// [[Rcpp::export]]
double timeDepBranch(Rcpp::NumericVector observations, std::string file, Rcpp::NumericVector initial, Rcpp::List transitions, Rcpp::List stops, bool silence, SEXP seed = R_NilValue){

	double seedcpp;
	if(Rf_isNull(seed)){
		seedcpp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	} else {
		seedcpp = Rf_asReal(seed);
	}
    gsl_rng_set(rng, seedcpp);
	silent = silence;


	#ifdef _WIN32
	HINSTANCE lib_handle;
	#else
	std::vector<void*> handles;
	#endif


	if(!silent) std::cout << "Starting process... " << std::endl;
	int nTrans = transitions.length();
	int nStops = stops.length();

	if(!silent) std::cout << "Initialization system..." << std::endl;
	// Initial population sizes
	std::vector<long int> init(initial.begin(), initial.end());

	// Initialize system
	System sys(init);

	// Add transitions
	if(!silent) std::cout << "Adding transitions..." << std::endl;

	// Iterate over transitions list
	for(int i = 0; i < nTrans; i++){

		// Get the ith transition
		Rcpp::List list_i = Rcpp::as<Rcpp::List>(transitions[i]);

		// Get popuation, is_random, probability
		int population = list_i[0];
		bool is_random = list_i[1];
		Rcpp::List rateList = Rcpp::as<Rcpp::List>(list_i[2]);
		Rate* r = nullptr;
		if(rateList.containsElementNamed("type")){
			if(Rcpp::as<int>(rateList["type"]) == 0)
			{

				if(!silent) std::cout << "Constant rate found!" << std::endl;
				r = new ConstantRate(Rcpp::as<double>(rateList[1]));

			} else if (Rcpp::as<int>(rateList["type"]) == 1)
			{

				Rcpp::NumericVector params = Rcpp::as<Rcpp::NumericVector>(rateList["params"]);
				if(!silent) std::cout << "Linear rate found!" << std::endl;
				r = new LinearRate(params[0], params[1]);

			} else if (Rcpp::as<int>(rateList["type"]) == 2){

				Rcpp::NumericVector params = Rcpp::as<Rcpp::NumericVector>(rateList["params"]);
				if(!silent) std::cout << "Switch rate found!" << std::endl;
				r = new SwitchRate(params[0], params[1], params[2]);

			} else if (Rcpp::as<int>(rateList["type"]) == 3){

				Rcpp::StringVector params = Rcpp::as<Rcpp::StringVector>(rateList["params"]);
				if(!silent) std::cout << "Custom rate found!" << std::endl;
				double (*rate)(double, void*);

				// Name for plugin library location
				const char* plugin_location = CHAR(Rf_asChar(params[0]));
				std::cout << params[0] << std::endl;
				std::cout << params[1] << std::endl;
				//Rcpp::Rcout << params[1] << std::endl;
				// Load plugin library
				#ifdef _WIN32
					lib_handle = LoadLibrary(plugin_location);
					if(!lib_handle)
					{
						Rcpp::stop("invalid file name for custom dll");
					}
				rate = (double (*)(double, void*))GetProcAddress(lib_handle, params[1]);
				#else
					void* hand = dlopen(plugin_location, RTLD_NOW);
					handles.push_back(hand);
					if(!hand)
					{
						Rcpp::stop("invalid file name custom dll");
					}
					rate = (double (*)(double, void*))dlsym(hand, params[1]);
				#endif

				r = new Rate(rate);
			}
			else{
				Rcpp::stop("invalid rate selection");
			}
		} else {
			r = new ConstantRate(Rcpp::as<double>(rateList[0]));
		}


		Update u;
		// If there is a random component, get which indices of the vector will be generated randomly
		if(is_random){
			// Get the offspring distribution of the Transition
			Rcpp::NumericVector oVec = Rcpp::as<Rcpp::NumericVector>(list_i[3]);
			std::vector<double> offspringVec (oVec.begin(), oVec.end());

			// Get the distribution
			std::string dist = list_i[4];

			// Get the offspring distribution of the Transition
			Rcpp::NumericVector p = Rcpp::as<Rcpp::NumericVector>(list_i[5]);
			std::vector<double> offspringParams (p.begin(), p.end());

			// Add this transition
			u = Update(is_random, offspringVec, dist, offspringParams);
		} else{
			// Get the fixed portion of the Transition
			Rcpp::NumericVector fix = Rcpp::as<Rcpp::NumericVector>(list_i[3]);
			std::vector<int> fixed (fix.begin(), fix.end());

			// add update
			u = Update(fixed);
		}

		sys.addUpdate(r, population, u);
	}

	// Add StopCriteria
	if(!silent) std::cout << "Adding stopping criteria..." << std::endl;

	// Iterate over StopCriterionList list
	for(int i = 0; i < nStops; i++){

		// Get the ith SC
		Rcpp::List list_i = Rcpp::as<Rcpp::List>(stops[i]);

		// Get the indices to combine, inequality for comparison, and value to compare against
		Rcpp::NumericVector ind = Rcpp::as<Rcpp::NumericVector>(list_i[0]);
		std::vector<int> indices (ind.begin(), ind.end());

		std::string ineq = list_i[1];
		double value = list_i[2];

		sys.addStop(StopCriterion(indices, ineq, value));
	}

	//sys.print();

	// Observation times
	std::vector<double> obsTimes(observations.begin(), observations.end());

	// Simulate
	if(!silent) std::cout << "Simulating..." << std::endl;
	sys.simulate_timedep(obsTimes, file);
	if(!silent) std::cout << "Ending process..." << std::endl;

	#ifdef _WIN32
		FreeLibrary(lib_handle);
	#else
		for(auto i = handles.begin(); i != handles.end(); i++){
			dlclose(*i);
			std::cout << "closing" << std::endl;
		}
	#endif

	return 0.0;
}
