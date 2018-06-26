/*
 * =====================================================================================
 *
 *       Filename:  State.cpp
 *
 *    Description:  Class representing system state
 *
 *        Version:  1.0
 *        Created:  06/13/20178 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University
 *
 * =====================================================================================
 */


#include "State.h"
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

State::State(){

}

State::State(std::vector<int> s){
	vec = s;

}

State::~State(){
	for(size_t i = 0; i < pops.size(); i++){
		delete pops[i];
	}
}

void State::print(){
	for(size_t i = 0; i < vec.size(); i++){
		std::cout << vec[i] << "\t";
	}
	std::cout << std::endl;
}

void State::toFile(int time, std::string file){
	std::ofstream of;

    of.open(file, std::fstream::in | std::fstream::out | std::fstream::app);
	of << time << "," << vec[0];
	for(size_t i = 1; i < vec.size(); i++){
		of << "," << vec[i];
	}

	of << std::endl;
}

void State::updateState(std::vector<int> update){
	if(update.size() != vec.size()){
		throw std::invalid_argument("States do not have same dimension");
	}

	for(size_t i = 0; i < update.size(); i++){
		vec[i] += update[i];
		if(vec[i] < 0)
			vec[i] = 0;
	}
}

void State::addPopulation(Population* p){
	pops.push_back(p);
}

Population* State::getPop(int i){
	return pops[i];
}

double State::getNextTime(){
	double overall = 0.0;
	for(size_t i = 0; i < pops.size(); i++){
		overall += vec[i] * pops[i]->getRate();
	}
	return(gsl_ran_exponential(rng, 1 / overall));
}

int State::choosePop(){
	std::vector<double> rates;
	for(size_t i  = 0; i < pops.size(); i++){
		rates.push_back(vec[i] * pops[i]->getRate());
		//std::cout << "Vec[i]: " << vec[i] << " Rate " << i << ": " << rates[i] << std::endl;
	}

	int choice = choose(rates);
	return choice;
}

void State::simulate(){
	std::cout << "In simulate()" << std::endl;
	double time = 0.0;
	for(int i = 0; i < 1000; i ++){
		double toNext = getNextTime();
		int pop = choosePop();
		std::vector<int> update = pops[pop]->getUpdate(gsl_rng_uniform(rng));
		updateState(update);
		std::cout << "Time: " << time + toNext << std::endl;
		print();
		time += toNext;
	}
}

void State::simulate(int numTime, std::string file){

	bool verbose = true;

	// Set up vector of observation times
    std::vector<int> obsTimes;
    for (int k = 0; k < numTime + 1; k++)
    {
        obsTimes.push_back(k);
    }
    //std::cout << "numTime: " << numTime << std::endl;

    // Set variables to keep track of our current time and which observation time comes next
    double curTime = 0;
    int curObsIndex = 0;

	int obsMod = pow(10, round(log10(numTime)-1));

    // Display some stuff - debugging
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
		/*
        if(checkZero())
        {
            std::cout << "A population has size 0, ending simulation" << std::endl;
            std::cout << "End time: " << curTime << std::endl;
            break;
        }
		*/

        // Get the next event time
        double timeToNext = getNextTime();
		//std::cout << "Time to next event: " << curTime + timeToNext << std::endl;


        // If our next event time is later than observation times,
        // Make our observations
        while((curTime + timeToNext > obsTimes[curObsIndex]))// & (curTime + timeToNext <= obsTimes[numTime]))
        {
			//print();
			toFile(obsTimes[curObsIndex], file);

			if(verbose &&  obsTimes[curObsIndex] % obsMod == 0)
				std::cout << "Time " << obsTimes[curObsIndex] << " of " << numTime << std::endl;
      curObsIndex++;
            //std::cout << "Sucessfully made an observation" << std::endl;
			if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;
        }
		if((unsigned)curObsIndex >= obsTimes.size()-1)
				break;

        // Update our state
		//std::cout << "Updating state..." << std::endl;
        int pop = choosePop();
		//std::cout << "Pop: " << pop << std::endl;
		std::vector<int> update = pops[pop]->getUpdate(gsl_rng_uniform(rng));
		/*for(size_t i = 0; i < update.size(); i++){
			std::cout << update[i] << " ";
		}
		std::cout << std::endl;
		*/
		updateState(update);
		//std::cout << "Finished updating state..." << std::endl;
		

        // Increase our current time and get the next Event Time
        curTime = curTime + timeToNext;
    }
	//std::cout << "curObsIndex: " << curObsIndex << std::endl;
	std::cout << "End Simulation Time: " << obsTimes[curObsIndex] << std::endl;
}
