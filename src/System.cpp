/*
 * =====================================================================================
 *
 *       Filename:  System.cpp
 *
 *    Description:  Class representing system
 *
 *        Version:  1.0
 *        Created:  06/13/2018 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */


#include "System.h"
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

System::System(){

}

System::System(std::vector<int> s){
	state = s;

}

System::~System(){
}

void System::print(){
	for(size_t i = 0; i < state.size(); i++){
		std::cout << state[i] << "\t";
	}
	std::cout << std::endl;
}

void System::toFile(int time, std::string file){
	std::ofstream of;

    of.open(file, std::fstream::in | std::fstream::out | std::fstream::app);
	of << time << "," << state[0];
	for(size_t i = 1; i < state.size(); i++){
		of << "," << state[i];
	}

	of << std::endl;
}

void System::updateSystem(std::vector<int> update){
	if(update.size() != state.size()){
		throw std::invalid_argument("Systems do not have same dimension");
	}

	for(size_t i = 0; i < update.size(); i++){
		state[i] += update[i];
		if(state[i] < 0)
			state[i] = 0;
	}
}

void System::addUpdate(double r, int f, Update u){
	rates.push_back(r);
	from.push_back(f);
	updates.push_back(u);
}

double System::getNextTime(std::vector<double>& o_rates){
	//double overall = 0.0;
	for(size_t i = 0; i < rates.size(); i++){
		//overall += rates[i] * state[i];
		o_rates[i] = rates[i] * state[from[i]];
	}
	return(gsl_ran_exponential(rng, 1 / std::accumulate(o_rates.begin(), o_rates.end(), 0.0)));
}

void System::simulate(int numTime, std::string file){
	bool verbose = true;

	std::vector<double> o_rates;
	for(size_t i = 0; i < rates.size(); i++){
		o_rates.push_back(0.0);
	}

	// Set up vector of observation times
    std::vector<int> obsTimes;
    for (int k = 0; k < numTime + 1; k++)
    {
        obsTimes.push_back(k);
    }


    // Set variables to keep track of our current time and which observation time comes next
    double curTime = 0;
    int curObsIndex = 0;

	int obsMod = pow(10, round(log10(numTime)-1));

    // Display some stuff  if verbose
	if(verbose){
		std::cout << "numTime: " << numTime << std::endl;
		std::cout << "Simulation Start Time: " << curTime << std::endl;
		std::cout << "Simulation End Time: " << obsTimes[numTime] << std::endl;
		std::cout << "obsTimes.size(): " << obsTimes.size() << std::endl;
	}
	
	//int wait = 0;

    // Run until our currentTime is greater than our largest Observation time
    while(curTime <= obsTimes[numTime])
    {
        Rcpp::checkUserInterrupt();
		
		bool zero = true;
		for(size_t i = 0; i < state.size(); i++){
			if(state[i] > 0)
				zero = false;
		}
		
		if(zero){
			std::cout << "All populations have gone extinct.  Exiting simulation..." << std::endl;
			break;
		}

        // Get the next event time
        double timeToNext = getNextTime(o_rates);
		/*for(size_t i = 0; i < o_rates.size(); i++){
			std::cout << "o_rates[" << i << "]: " << o_rates[i] << std::endl;
		}
		
		if(wait == 0){
			Rcpp::Environment base = Rcpp::Environment("package:base");
			Rcpp::Function readline = base["readline"];
			Rcpp::Function as_numeric = base["as.numeric"];
			wait = Rcpp::as<int>(as_numeric(readline("> ")));
		}
		std::cout << "Time to next event: " << curTime + timeToNext << std::endl;
		*/

        // If our next event time is later than observation times,
        // Make our observations
        while((curTime + timeToNext > obsTimes[curObsIndex]))// & (curTime + timeToNext <= obsTimes[numTime]))
        {
			// print out current state vector
			toFile(obsTimes[curObsIndex], file);

			if(verbose &&  obsTimes[curObsIndex] % obsMod == 0)
				std::cout << "Time " << obsTimes[curObsIndex] << " of " << numTime << std::endl;
      curObsIndex++;
			if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;
        }
		if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;

        // Update our System
		//std::cout << "Updating System..." << std::endl;
        int index = choose(o_rates);

		//std::cout << "Choice: " << index << std::endl;

		std::vector<int> update = updates[index].get();
		update[index] = update[index] - 1;
		updateSystem(update);
		//std::cout << "Finished updating System..." << std::endl;


        // Increase our current time and get the next Event Time
        curTime = curTime + timeToNext;
    }
	//std::cout << "curObsIndex: " << curObsIndex << std::endl;
	std::cout << "End Simulation Time: " << obsTimes[curObsIndex] << std::endl;
}
