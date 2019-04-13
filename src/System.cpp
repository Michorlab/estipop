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
#include <sstream>
#include <iomanip>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;
extern bool silent;

void out(std::string s){
	//std::cout << s << std::endl;
}

template <typename T>
std::string to_string_wp(const T a_value, const int n = 20)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

System::System(){

}

System::System(std::vector<int> s){
	state = s;
	tot_rate_homog = 0.0;
}

System::~System(){
	for(size_t i = 0; i < rates2.size(); i++){
		delete rates2[i];
	}
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

void System::addUpdate(Rate* r, int f, Update u){
	rates2.push_back(r);
	from.push_back(f);
	updates.push_back(u);

	tot_rate_homog += r->rate_homog;
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
    while(curTime <= obsTimes[obsTimes.size()-1])
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
	std::cout.precision(20);
	double tot_rate;
	double rand_next_time = 0.0;

	tot_rate_homog = 0.0;

	for(size_t i = 0; i < rates2.size(); i++){
		tot_rate_homog += rates2[i]->rate_homog * state[from[i]];
	}

	out("In nextTime2");

	while(true){
		Rcpp::checkUserInterrupt();

		tot_rate = 0;

		out("tot_rate_homog: " + to_string_wp(tot_rate_homog));
		rand_next_time += gsl_ran_exponential(rng, 1 / tot_rate_homog);


		out("starting loop");
		for(size_t i = 0; i < rates2.size(); i++){
			tot_rate += (*rates2[i])(curTime + rand_next_time) * state[from[i]];
		}
		out("ending loop");

		double u_thin = gsl_ran_flat(rng, 0, 1);
		double beta_ratio = tot_rate / tot_rate_homog;
		//std::cout << beta_ratio << "\n";

		if(u_thin <= beta_ratio)
		{
		  break;
		}
	}

	out("Returning from nextTime2 with: " + to_string_wp(rand_next_time));
	if(rand_next_time <= 0){
		out("I'VE RETURNED A ZERO TIME TO NEXT EVENT!!!");
	}

	return(rand_next_time);
}

int System::getNextEvent2(double curTime, double timeToNext){
	std::vector<double> cumulativeHazards;

	out("in nextEvent2");


	for(size_t i = 0; i < rates2.size(); i++){
		cumulativeHazards.push_back((*rates2[i])(curTime, curTime + timeToNext) * state[from[i]]);
	}

	out("returning from nextEvetn2");
	return choose(cumulativeHazards);
}

void System::simulate2(int numTime, std::string file){
	bool verbose = false;
	std::cout.precision(60);

	std::vector<double> o_rates;
	for(size_t i = 0; i < rates.size(); i++){
		o_rates.push_back(0.0);
	}

	/*
	// Set up vector of observation times
    std::vector<int> obsTimes;
    for (int k = 0; k < numTime + 1; k++)
    {
        obsTimes.push_back(k);
    }
	*/
	std::vector<double> obsTimes;
    for (double k = 0; k <= numTime; k+=1.00)
    {
        obsTimes.push_back(k);
    }
	numTime = obsTimes.size()-1;



    // Set variables to keep track of our current time and which observation time comes next
    double curTime = 0;
    int curObsIndex = 0;

	int obsMod = std::max(1, (int)pow(10, round(log10(numTime)-1)));

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

		out("Got my next time: " + to_string_wp(timeToNext));

        // If our next event time is later than observation times,
        // Make our observations
        while((curTime + timeToNext > obsTimes[curObsIndex]) & curObsIndex <= obsTimes.size()-1)// & (curTime + timeToNext <= obsTimes[numTime]))
        {
			out("Time to next: " + to_string_wp(curTime + timeToNext) + " curObsIndex: " + to_string_wp(curObsIndex));
			out("In the observation loop");
			// print out current state vector
			toFile(obsTimes[curObsIndex], file);

			out("Checking this");
			if((int)obsTimes[curObsIndex] % obsMod == 0){
				std::cout << "Time " << (int)obsTimes[curObsIndex] << " of " << numTime << std::endl;
			}

			curObsIndex++;

			out("Checking here");
			/*if((unsigned)curObsIndex >= obsTimes.size()-1){
				out("I'm quitting in here: " + to_string_wp(curObsIndex) + ", " + to_string_wp(obsTimes.size()-1));
				break;
			}*/
        }

		if((unsigned)curObsIndex > obsTimes.size()-1){
			out("I'm quitting in here: " + to_string_wp(curObsIndex) + ", " + to_string_wp(obsTimes.size()-1));
				break;
		}


        // Update our System
		out("headed from sim2 to getEvent2");
        int index = getNextEvent2(curTime, timeToNext);


		std::vector<int> update = updates[index].get();
		update[from[index]] = update[from[index]] - 1;
		updateSystem(update);

        // Increase our current time and get the next Event Time
        curTime = curTime + timeToNext;
		out("Curtime: " + to_string_wp(curTime));

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
	//toFile(obsTimes[obsTimes.size() - 1], file);
	std::cout << "End Simulation Time: " << obsTimes[curObsIndex-1] << std::endl;
}
