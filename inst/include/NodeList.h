/*
 * =====================================================================================
 *
 *       Filename:  NodeList.h
 *
 *    Description:  Class representing a list of cell states and count
 *
 *        Version:  1.0
 *        Created:  03/01/2017 14:57:27
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

#include "Node.h"


class NodeList {
public:
	// Members

	Node* root;
	Node* cur;

	int total;
	int left;
	int right;

	double totalFitness;

	int numMutated;

	// Constructure / Destructor
	NodeList();
	~NodeList();

	// Manipulations, including methods for holding DiffTriangle information
	void insert(long int b, std::string m, double f, int c);
	void insert(long int b, std::string m, double f, int tri, int l, int c);
	void remove(long int b, std::string m, double f, int c);
	void remove(long int b, std::string m, double f, int tri, int l, int c);

	// Print methods
	void traversePrint();
	void printCur();

	// Methods to keep tree balanced by count
	void balanceLeft();
	void balanceRight();

	// Indexing methods
	Node* getAt(int index);
	void removeAt(int index);
	void increaseAt(int index);

	// Fitness indexing
	Node* getAtFitness(double f);

	// Methods using cur, used immediately after cur has been set to save computation time
	Node* getCur();
	void setCurIndex(int index);
	void setCurFitness(double f);
	void increaseAtCur();
	void removeAtCur();

	// Outputs
	double labelled();
	double sdi();
	std::string diversity();
	double shannon(std::map<long int, int> c_map);
	double simpson(std::map<long int, int> c_map);

	std::map<long int, int> count_map_codes();
	std::map<std::string, int> count_map_mutation();
	int numTypes();

	// Free memory
	void deleteList();

	// Two methods to write to file
	void writeToFile(std::ofstream& of, int time);
	void writeToFile2(std::ofstream& of, int time);

	// Fast method to build a DiffTriangle's NodeList structure
	std::vector<std::vector<Node*>>makeTriangle(int nlevels, int mfac);
};

