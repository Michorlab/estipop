/*
 * =====================================================================================
 *
 *       Filename:  State.h
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

#pragma once
#include <string>
#include <iostream>
#include <ostream>
#include <sstream>
#include <vector>
#include <map>

#include "Population.h"

class State {
public:
	// Members
	std::vector<int> vec;
	std::vector<Population*> pops;
	
	// Constructors
	State();
	State(std::vector<int> s);
	~State();
	
	// Methods
	void print();
	void toFile(int time, std::string file);
	void updateState(std::vector<int> update);
	
	void addPopulation(Population* p);
	Population* getPop(int i);
	
	double getNextTime();
	int choosePop();
	
	void simulate();
	void simulate(int numTime, std::string file);
};

