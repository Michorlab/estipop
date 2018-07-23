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

class System {
public:
	// Members
	std::vector<int> state;
	
	std::vector<double> rates;
	std::vector<int> from;
	std::vector<Update> updates;
	
	// Constructors
	System();
	System(std::vector<int> s);
	~System();
	
	// Methods
	void print();
	void toFile(int time, std::string file);
	void updateSystem(std::vector<int> update);
	
	void addUpdate(double r, int f, Update u);
	
	double getNextTime(std::vector<double>& o_rates);
	
	void simulate(int numTime, std::string file);
};

