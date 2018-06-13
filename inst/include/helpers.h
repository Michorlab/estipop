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

// Helper methods - CellPopulationCode
long int cellByFitness(std::vector<double> fitness);
std::vector<double> normalize(std::vector<double> input);
int choose(std::vector<double> input);
std::map<long int, int> count_map(std::vector<long int> input);
std::map<int, int> count_map(std::vector<int> input);
void printVec(std::vector<long int> v);
void printVec(std::vector<double> v);

// Helper methods - DiffTriCode
int beginIndex(int mfac, int level);
std::vector<std::vector<long int> > splitDoubleVector(std::vector<long int> v, std::vector<double> probs);



