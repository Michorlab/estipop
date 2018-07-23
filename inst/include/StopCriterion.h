/*
 * =====================================================================================
 *
 *       Filename:  StopCriterion.h
 *
 *    Description:  Class representing a system StopCriterion
 *
 *        Version:  1.0
 *        Created:  06/26/20178 14:57:27
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
//#include <map>



class StopCriterion {
public:
	std::vector<int> indices;
	std::string inequality;
	double value;
	
	StopCriterion();
	StopCriterion(std::vector<int> ind, std::string ineq, double v);
	~StopCriterion();
	
	bool check(std::vector<int> state);
};

