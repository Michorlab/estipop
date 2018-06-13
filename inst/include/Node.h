/*
 * =====================================================================================
 *
 *       Filename:  Node.h
 *
 *    Description:  Class representing a single cell state and count
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
#include <vector>
#include <map>

class NodeList;

class Node {
public:
	// Members
	
	Node* next;
	Node* prev;
	NodeList* list;
	
	// Information
	long int barcode;
	std::string mutation;
	double fitness;
	
	// Triangle Information
	int triangle;
	int level;
	
	int count;
	
	// Counstructors
	Node();
	Node(long int b, std::string m, double f);
	Node(long int b, std::string m, double f, int c);
	Node(long int b, std::string m, double f, int t, int l);
	Node(long int b, std::string m, double f, int t, int l, int c);
	
	// Used to associate a Node to its NodeList
	void setList(NodeList* l);
	
	// Remove a Node from a list
	void removeFromList();
	
	// Output
	void printNode();
};

