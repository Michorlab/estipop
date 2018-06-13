/*
 * =====================================================================================
 *
 *       Filename:  Node.cpp
 *
 *    Description:  Class representing a single cell state and count
 *
 *        Version:  1.0
 *        Created:  03/01/2017 14:57:27
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jeremy Ferlic (), jferlic@g.harvard.edu
 *   Organization:  Harvard University,
 *
 * =====================================================================================
 */

#include "Node.h"
#include "NodeList.h"
#include <iostream>
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

Node::Node(){
	next = NULL;
	prev = NULL;
}

Node::Node(long int b, std::string m, double f){
	next = NULL;
	prev = NULL;
	
	barcode = b;
	mutation = m;
	fitness = f;
	
	triangle = -1;
	level = -1;
}

Node::Node(long int b, std::string m, double f, int c){
	next = NULL;
	prev = NULL;
	
	barcode = b;
	mutation = m;
	fitness = f;
	count = c;
	
	triangle = -1;
	level = -1;
}

Node::Node(long int b, std::string m, double f, int t, int l){
	next = NULL;
	prev = NULL;
	
	barcode = b;
	mutation = m;
	fitness = f;
	
	triangle = t;
	level = l;
}

Node::Node(long int b, std::string m, double f, int t, int l, int c){
	next = NULL;
	prev = NULL;
	
	barcode = b;
	mutation = m;
	fitness = f;
	
	triangle = t;
	level = l;
	
	count = c;
}

void Node::setList(NodeList* l){
	list = l;
}

void Node::removeFromList(){
	if(list->root == this){
		list->root = this->next;
	}
	list->total -= this->count;
	list->totalFitness -= this->count * this->fitness;
	
	if(this->mutation != "0")
		list->numMutated -= this->count;
	
	if(this->prev != NULL)
		this->prev->next = this->next;
	
	if(this->next != NULL)
		this->next->prev = this->prev;
	
	delete this;
}

void Node::printNode(){
	std::cout << "Code: " << barcode << " Mutation: " << mutation << " Fitness: " << fitness 
	<< " Triangle: " << triangle << " Level: " << level << " Count: " << count << std::endl;
}