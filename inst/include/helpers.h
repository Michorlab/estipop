/*
 * =====================================================================================
 *
 *       Filename:  helpers.h
 *
 *    Description:  Helper Functions
 *
 *        Version:  1.0
 *        Created:  05/01/2018 14:57:27
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
#include <vector>
#include <map>
#include <gsl/gsl_math.h>

// Helper methods - CellPopulationCode
std::vector<double> normalize(std::vector<double> input);
int choose(std::vector<double> input);

// Rate functions
double maximizeFunc(gsl_function rate_function, double start_time, double end_time, int bins);

void maximizePiecewise(gsl_function rate_function, double start_time, double end_time, int bins, std::vector<double>& vec, double buffer);




