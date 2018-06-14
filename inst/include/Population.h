/*
 * =====================================================================================
 *
 *       Filename:  Population.h
 *
 *    Description:  Class representing a population
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



class Population {
public:
	// Members
	double rate;
	std::vector<double> probs;
	std::vector<double> norm_probs;
	std::vector<std::vector<int>> updates;
	
	// Constructors
	Population();
	Population(double r);
	~Population();
	
	// Methods
	double getRate();
	void setRate(double r);
	
	void addUpdate(double prob, std::vector<int> update);
	
	void normalizeProbs();
	
	std::vector<int> getUpdate(double prob);
	
	void printUpdates();
};

