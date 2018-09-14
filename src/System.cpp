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
#include <gsl/gsl_randist.h>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;
extern bool silent;

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

void System::toFile(double time, std::string file){
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

void System::addStop(StopCriterion c){
	stops.push_back(c);
}

double System::getNextTime(std::vector<double>& o_rates){
	for(size_t i = 0; i < rates.size(); i++){
		o_rates[i] = rates[i] * state[from[i]];
	}
	return(gsl_ran_exponential(rng, 1 / std::accumulate(o_rates.begin(), o_rates.end(), 0.0)));
}

void System::simulate(int numTime, std::string file){
	bool verbose = false;

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
	if(!silent){
		std::cout << "numTime: " << numTime << std::endl;
		std::cout << "Simulation Start Time: " << curTime << std::endl;
		std::cout << "Simulation End Time: " << obsTimes[numTime] << std::endl;
		std::cout << "obsTimes.size(): " << obsTimes.size() << std::endl;
	}

    // Run until our currentTime is greater than our largest Observation time
    while(curTime <= obsTimes[numTime])
    {
        Rcpp::checkUserInterrupt();

        // Get the next event time
        double timeToNext = getNextTime(o_rates);

        // If our next event time is later than observation times,
        // Make our observations
        while((curTime + timeToNext > obsTimes[curObsIndex]))// & (curTime + timeToNext <= obsTimes[numTime]))
        {
			// print out current state vector
			toFile(obsTimes[curObsIndex], file);

			if(verbose &&  obsTimes[curObsIndex] % obsMod == 0 && !silent)
				std::cout << "Time " << obsTimes[curObsIndex] << " of " << numTime << std::endl;
      curObsIndex++;
			if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;
        }
		if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;

        // Update our System
        int index = choose(o_rates);


		std::vector<int> update = updates[index].get();
		update[from[index]] = update[from[index]] - 1;
		updateSystem(update);

        // Increase our current time and get the next Event Time
        curTime = curTime + timeToNext;
		
		bool stop = false;

		for(size_t i = 0; i < stops.size(); i++){
			if(stops[i].check(state))
				stop = true;
		}

		if(stop){
			toFile(curTime, file);
			if(!silent)
				std::cout << "A stopping criterion has been met. Exiting simulation..." << std::endl;
			break;
		}

		bool zero = true;
		for(size_t i = 0; i < state.size(); i++){
			if(state[i] > 0)
				zero = false;
		}

		if(zero){
			toFile(curTime, file);
			if(!silent)
				std::cout << "All populations have gone extinct.  Exiting simulation..." << std::endl;
			break;
		}

    }
		if(!silent)
	std::cout << "End Simulation Time: " << obsTimes[curObsIndex] << std::endl;
}

double System::getNextTime2(double curTime){
	double tot_rate;
	double rand_next_time;
	
	while(true){
		Rcpp::checkUserInterrupt();
		
		tot_rate = 0;
		rand_next_time = gsl_ran_exponential(rng, 1 / tot_rate_homog);
		
		for(size_t i = 0; i < rates2.size(); i++){
			tot_rate += rates2[i].eval(curTime + rand_next_time) * state[from[i]];
		}
		
		double u_thin = gsl_ran_flat(rng, 0, 1);
		double beta_ratio = tot_rate / tot_rate_homog;
		//std::cout << beta_ratio << "\n";

		if(u_thin <= beta_ratio)
		{
		  break;
		}
	}

	return(rand_next_time);
}

int System::getNextEvent2(double curTime, double timeToNext){
	std::vector<double> cumulativeHazards;
	
	for(size_t i = 0; i < rates2.size(); i++){
		cumulativeHazards.push_back(rates2[1].integrateFunct(curTime, curTime + timeToNext));
	}
	
	return choose(cumulativeHazards);
}

void System::simulate2(int numTime, std::string file){
	bool verbose = false;

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

    // Run until our currentTime is greater than our largest Observation time
    while(curTime <= obsTimes[numTime])
    {
        Rcpp::checkUserInterrupt();

        // Get the next event time
        double timeToNext = getNextTime2(curTime);

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
        int index = getNextEvent2(curTime, timeToNext);


		std::vector<int> update = updates[index].get();
		update[from[index]] = update[from[index]] - 1;
		updateSystem(update);

        // Increase our current time and get the next Event Time
        curTime = curTime + timeToNext;
		
		bool stop = false;

		for(size_t i = 0; i < stops.size(); i++){
			if(stops[i].check(state))
				stop = true;
		}

		if(stop){
			toFile(curTime, file);
			std::cout << "A stopping criterion has been met. Exiting simulation..." << std::endl;
			break;
		}

		bool zero = true;
		for(size_t i = 0; i < state.size(); i++){
			if(state[i] > 0)
				zero = false;
		}

		if(zero){
			toFile(curTime, file);
			std::cout << "All populations have gone extinct.  Exiting simulation..." << std::endl;
			break;
		}

    }
	std::cout << "End Simulation Time: " << obsTimes[curObsIndex] << std::endl;
}
