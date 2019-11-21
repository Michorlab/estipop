/*
 * =====================================================================================
 *
 *       Filename:  System.h
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

#pragma once
#include <string>
#include <iostream>
#include <ostream>
//#include <sstream>
#include <vector>
//#include <map>

#include "Update.h"
#include "Rate.h"
#include "StopCriterion.h"

class System {
public:
	// Members
	double tot_rate_homog;

	std::vector<long int> state;

	std::vector<double> rates;
	std::vector<Rate*> rates2;
	std::vector<int> from;
	std::vector<Update> updates;

	std::vector<StopCriterion> stops;

	// Constructors
	System();
	System(std::vector<long int> s);
	~System();

	// Methods
	void print();
	void toFile(double time, std::string file);
	void updateSystem(std::vector<int> update);

	void addUpdate(double r, int f, Update u);

	void addUpdate(Rate* r, int f, Update u);

	void addStop(StopCriterion c);

	double getNextTime(std::vector<double>& o_rates);

	double getNextTime2(double curTime, double totTime);

	int getNextEvent2(double curTime, double timeToNext);

	void simulate(std::vector<double> obsTimes, std::string file);

	void simulate_timedep(std::vector<double> obsTimes, std::string file);
};
