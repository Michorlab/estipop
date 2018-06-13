/*
 * =====================================================================================
 *
 *       Filename:  NodeList.cpp
 *
 *    Description:  Class representing a list of cell states and count
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

#include "NodeList.h"

#include <iostream>
#include <fstream>
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

NodeList::NodeList(){
	root = NULL;
	cur = NULL;
	total = 0;
	left = 0;
	right = 0;
	totalFitness = 0;
	numMutated = 0;
}

NodeList::~NodeList(){
	deleteList();
}

void NodeList::insert(long int b, std::string m, double f, int c){
	if(m != "0"){
		numMutated += c;
	}
	if (root == NULL){
		cur = new Node(b, m, f, c);
		cur->next = NULL;
		cur->prev = NULL;
		cur->setList(this);
		root = cur;
	} else {
		cur = root;
		Node * t = root;
		left = 0;
		right = total - t->count;
		while(t != NULL){
			if(t->barcode == b && t->mutation == m && t->fitness == f){
				t->count += c;
				total += c;
				totalFitness += f*c;
				balanceLeft();
				return;
			}

			left += t->count;
			t = t->next;
			if(t != NULL){
				right -= t->count;
				cur = t;
			}
		}

		t = new Node(b, m, f, c);
		cur->next = t;
		t->prev = cur;
		cur = t;
		cur->setList(this);
		balanceLeft();
	}

	total += c;
	totalFitness += f*c;
}

void NodeList::insert(long int b, std::string m, double f, int tri, int l, int c){
	if(m != "0"){
		numMutated += c;
	}
	if (root == NULL){
		cur = new Node(b, m, f, tri, l, c);
		cur->next = NULL;
		cur->prev = NULL;
		cur->setList(this);
		root = cur;
	} else {
		cur = root;
		Node * t = root;
		left = 0;
		right = total - t->count;
		while(t != NULL){
			if(t->barcode == b && t->mutation == m && t->fitness == f && t->triangle == tri && t->level == l){
				t->count += c;
				total += c;
				totalFitness += f*c;
				balanceLeft();
				return;
			}

			left += t->count;
			t = t->next;
			if(t != NULL){
				right -= t->count;
				cur = t;
			}
		}

		t = new Node(b, m, f, tri, l, c);
		cur->next = t;
		t->prev = cur;
		cur = t;
		cur->setList(this);
		balanceLeft();
	}

	total += c;
	totalFitness += f*c;
}

void NodeList::remove(long int b, std::string m, double f, int c){
	if (m != "0"){
		numMutated -= c;
	}
	if (root == NULL){
		return;
	} else {
		cur = root;
		Node * t = root;
		left = 0;
		right = total - t->count;
		while(t != NULL){
			if(t->barcode == b && t->mutation == m && t->fitness == f){
				if(t->count - c <= 0){
					t->removeFromList();
					return;
				} else {
					t->count -= c;
					total -= c;
					totalFitness -= f*c;
					balanceRight();
					return;
				}
			}

			left += t->count;
			t = t->next;
			if(t != NULL){
				right -= t->count;
				cur = t;
			}
		}
	}
}

void NodeList::remove(long int b, std::string m, double f, int tri, int l, int c){
	if (m != "0"){
		numMutated -= c;
	}
	if (root == NULL){
		return;
	} else {
		cur = root;
		Node * t = root;
		left = 0;
		right = total - t->count;
		while(t != NULL){
			if(t->barcode == b && t->mutation == m && t->fitness == f && t->triangle == tri && t->level == l){
				if(t->count - c <= 0){
					t->removeFromList();
					return;
				} else {
					t->count -= c;
					total -= c;
					totalFitness -= f*c;
					balanceRight();
					return;
				}
			}

			left += t->count;
			t = t->next;
			if(t != NULL){
				right -= t->count;
				cur = t;
			}
		}
	}
}

void NodeList::balanceLeft(){
	Node* keep = cur;

	Node* t = cur;
	while(t->prev != NULL && t->prev->count < keep->count){
		t = t->prev;
		right += t->count;
		left -= t->count;
	}

	if(t != keep){
        // if keep is not end of list
		if(keep->next != NULL)
			keep->next->prev = keep->prev;
		keep->prev->next = keep->next;

        // if t is not the beginning of the list
		if(t->prev != NULL){
			t->prev->next = keep;
		} else {
			root = keep;
		}
		keep->prev = t->prev;
		keep->next = t;
		t->prev = keep;
	}
}

void NodeList::balanceRight(){
	Node* keep = cur;

	Node* t = cur;
	while(t->next != NULL && t->next->count > keep->count){
		t = t->next;
		right -= t->count;
		left += t->count;
	}

	if(t != keep){
        //change root if necessary
        if(root == keep){
            root = keep->next;
        }
        // stitch out keep
        if(keep->prev != NULL)
            keep->prev->next = keep->next;
        if(keep->next != NULL)
            keep->next->prev = keep->prev;

        // fix keeps prev
        keep->prev = t;
        keep->next = t->next;

        if(t->next != NULL)
            t->next->prev = keep;
        t->next = keep;
	}
}


Node* NodeList::getAt(int index){
	if(index >= total)
		return NULL;
	cur = root;
	left = 0;
	right = total - cur->count;
	int t = cur->count;
	while(t < index+1){
		cur = cur->next;
		t = t + cur->count;
	}
	//std::cout << "total: " << total << std::endl;
	return cur;
}

Node* NodeList::getCur(){
	return cur;
}

void NodeList::setCurIndex(int index){
	if(index >= total)
		throw "Invalid index!";
	cur = root;
	left = 0;
	right = total - cur->count;
	int t = cur->count;
	while(t < index+1){
		cur = cur->next;
		t = t + cur->count;
	}
	//std::cout << "total: " << total << std::endl;
	return;
}


void NodeList::setCurFitness(double f){
	if(f >= totalFitness){
		throw "Invalid fitness!";
	}
	cur = root;
	double curFit = cur->fitness * cur->count;

	while(curFit < f){
		cur = cur->next;
		curFit += (cur->count * cur->fitness);
	}

	return;
}

void NodeList::increaseAtCur(){
	Node* hold = cur;
	Node* t = cur;
	t->count += 1;
	total += 1;
	totalFitness += t->fitness;
	if(t->mutation != "0"){
		numMutated++;
	}
	balanceLeft();
	cur = hold;
}


void NodeList::removeAtCur(){
	Node* hold = cur;
	Node* t = cur;
	t->count -= 1;
	total -= 1;
	totalFitness -= t->fitness;
	cur = t;
	if(t->mutation != "0"){
		numMutated--;
	}
	balanceRight();
	cur = hold;
}

void NodeList::removeAt(int index){
	//std::cout << "In removeAt: " << index << std::endl;
	Node* t = getAt(index);

	if(t->count - 1 <= 0){
		t->removeFromList();
		return;
	} else {
		t->count -= 1;
		total -= 1;
		totalFitness -= t->fitness;
		cur = t;
		if(t->mutation != "0"){
			numMutated--;
		}
		/*if(t->count == 0){
			t->prev->next = t->next;
			t->next->prev = t->prev;
			if(root == t)
				root = t->next;
			delete t;
		}
		cur = root;
		*/
		//std::cout << "about to head to balanceRight()" << std::endl;
		balanceRight();
	}
}

void NodeList::increaseAt(int index){
	Node* t = getAt(index);
	t->count += 1;
	total += 1;
	totalFitness += t->fitness;
	cur = t;
	if(t->mutation != "0"){
		numMutated++;
	}
	balanceLeft();
}

Node* NodeList::getAtFitness(double f){
	if(f >= totalFitness){
		return NULL;
	}
	cur = root;
	double curFit = cur->fitness * cur->count;

	while(curFit < f){
		cur = cur->next;
		curFit += (cur->count * cur->fitness);
	}

	return cur;
}

void NodeList::traversePrint(){
	for(Node* t = root; t != NULL; t = t->next){
		std::cout << "Code: " << t->barcode << " Mutation: " << t->mutation << " Fitness: " << t->fitness
	<< " Triangle: " << t->triangle << " Level: " << t->level << " Count: " << t->count << std::endl;
	}
	std::cout << "Total: " << total << std::endl;
	std::cout << "Types: " << numTypes() << std::endl;
	std::cout << "Total Fitness: " << totalFitness << std::endl;
	std::cout << "Number Mutated: " << numMutated << std::endl;
}

void NodeList::printCur(){
	std::cout << "Left: " << left << " Current: " << cur->count << " Right: " << right << " Total: " << total << std::endl;
}

double NodeList::labelled(){
	int labelled = 0;
	int total = 0;
	for(Node* t = root; t != NULL; t = t->next){
		if(t->barcode >= 0){
			labelled += t->count;
		}

		total += t->count;
	}

	return (double) labelled / (double) total;
}

double NodeList::sdi(){
	// Total number of codes N and counter map of codes
	long int N = total;
	std::map<long int, int> c_map = count_map_codes();

	// Shannon Diversity simply a sum of terms
	double ret = 0.0;
	std::map<long int, int>::iterator it;
	for ( it = c_map.begin(); it != c_map.end(); it++ )
	{
		double c = (double)it->second;
		ret = ret + (c / N * log(c / N));
	}
	return -ret;
}

std::string NodeList::diversity(){
	
	std::map<long int, int> c_map = count_map_codes();
	std::ostringstream o;
    o << c_map.size() << "," << shannon(c_map) << "," << simpson(c_map);
    return o.str();
}

double NodeList::shannon(std::map<long int, int> c_map){
	// Total number of codes N and counter map of codes
	long int N = total;
	
	// Shannon Diversity simply a sum of terms
	double ret = 0.0;
	std::map<long int, int>::iterator it;
	for ( it = c_map.begin(); it != c_map.end(); it++ )
	{
		double c = (double)it->second;
		ret = ret + (c / N * log(c / N));
	}
	return -ret;
}

double NodeList::simpson(std::map<long int, int> c_map){
	// Total number of codes N and counter map of codes
	long int N = total;

	// Shannon Diversity simply a sum of terms
	double ret = 0.0;
	std::map<long int, int>::iterator it;
	for ( it = c_map.begin(); it != c_map.end(); it++ )
	{
		double c = (double)it->second;
		ret = ret + ((c / N) * (c / N));
	}
	return 1/ret;
}

std::map<long int, int> NodeList::count_map_codes(){
	std::map<long int, int> m;
	for(Node* t = root; t != NULL; t = t->next){
		long int code = t->barcode;

		std::map<long int,int>::iterator it = m.find(code);
		if (it != m.end()){
			m[code] += t->count;
		} else {
			if(t->count > 0){
				m.insert({code, t->count});
			}
		}
	}

	return m;
}

std::map<std::string, int> NodeList::count_map_mutation(){
	std::map<std::string, int> m;
	for(Node* t = root; t != NULL; t = t->next){
		std::string mut = t->mutation;

		std::map<std::string,int>::iterator it = m.find(mut);
		if (it != m.end()){
			m[mut] += t->count;
		} else {
			if(t->count > 0){
				m.insert({mut, t->count});
			}
		}
	}

	return m;
}

int NodeList::numTypes(){
	int ret = 0;
	for(Node* t = root; t != NULL; t = t->next){
		if(t->count > 0)
			ret++;
	}

	return ret;
}

void NodeList::deleteList(){
	total = 0;
	left = 0;
	right = 0;
	Node* pnode;
	while (root != NULL)
	{
		pnode = root;
		root = root->next;
		delete pnode;
	}
}

void NodeList::writeToFile(std::ofstream& of, int time){
	for(Node* t = root; t != NULL; t = t->next){
		of << time << "," << t->barcode << "," << t->mutation << "," << t->fitness << "," << t->count << std::endl;
	}
	//of << "Total: " << total << std::endl;
	//of << "Types: " << numTypes() << std::endl;
	//of << "Total Fitness: " << totalFitness << std::endl;
	//of << "Number Mutated: " << numMutated << std::endl;
}

void NodeList::writeToFile2(std::ofstream& of, int time){
	NodeList temp;
	for(Node* t = root; t != NULL; t = t->next){
		temp.insert(t->barcode, t->mutation, t->fitness, t->count);
		//of << "C: " << t->barcode << " M: " << t->mutation << " F: " << t->fitness << " Count: " << t->count << std::endl;
	}
	temp.writeToFile(of, time);
	/*
	of << "Total: " << total << std::endl;
	of << "Types: " << numTypes() << std::endl;
	of << "Total Fitness: " << totalFitness << std::endl;
	of << "Number Mutated: " << numMutated << std::endl;
	*/
}

std::vector<std::vector<Node*>> NodeList::makeTriangle(int nlevels, int mfac){
  std::vector<std::vector<Node*>> ret;
	for(int i = 0; i < mfac; i++){
		std::vector<Node*> t;
		for(int j = 0; j < nlevels; j++){
			t.push_back(NULL);
		}
		ret.push_back(t);
	}

	if(total > 0)
	  return ret;

	Node* p;

    for(int j = nlevels; j > 0; j--)
    {
		for(int i = 0; i < mfac; i++)
		{
			if(i == 0 && j == nlevels)
			{
				Node* t = new Node();
				t->barcode = -1;
				t->mutation = "0";
				t->fitness = 1;
				t->triangle = i;
				t->level = j;
				t->count = std::pow(2, j-1);
				t->list = this;
				root = t;
				p = t;

				ret[i][j-1] = t;
			} else {
				Node* t = new Node();
				t->barcode = -1;
				t->mutation = "0";
				t->fitness = 1;
				t->triangle = i;
				t->level = j;
				t->count = std::pow(2, j-1);
				t->list = this;
				t->prev = p;
				p->next = t;
				p = t;

				ret[i][j-1] = t;
			}
        }
    }

	total = mfac * (std::pow(2, nlevels) - 1);
	cur = root;
	left = 0;
	right = total - std::pow(2, nlevels-1);
	totalFitness = total;

	return ret;
}

