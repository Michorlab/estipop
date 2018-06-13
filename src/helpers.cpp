/*
 * =====================================================================================
 *
 *       Filename:  helpers.cpp
 *
 *    Description:  implementing various helper methods
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
#include <gsl/gsl_randist.h>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;

// Given a list of fitness values, choose one according to the fitness values
long int cellByFitness(std::vector<double> fitness)
{
    gsl_ran_discrete_t * F;
    F = gsl_ran_discrete_preproc(fitness.size(), &fitness[0]);

    double draw = gsl_ran_discrete(rng, F);

    gsl_ran_discrete_free(F);
    return draw;
}

// Helper method to take vector, normalize to length 1, and then return the cumulative vector
std::vector<double> normalize(std::vector<double> input)
{
    double s = std::accumulate(input.begin(), input.end(), 0.0);
    //std::cout << "s: " << s << std::endl;
    std::transform(input.begin(), input.end(), input.begin(),
                   std::bind1st(std::multiplies<double>(),1/s));
    std::vector<double> sum(input.size());
    std::partial_sum(input.begin(), input.end(), sum.begin());
    return(sum);
}

// Helper method
// Input: list of 'n' doubles
// Output: a choice from 1 to n according to probability from input
int choose(std::vector<double> input)
{
    std::vector<double> c_norm = normalize(input);
    double r = gsl_rng_uniform(rng);
  
    // Choose which input
    for(size_t i = 0; i < input.size(); i++)
    {
        if(r <= c_norm[i])
        {
            return i;
        }

    }
    return -1;
}

// Helper method - return a count map of instances of a barcode in a cell population
std::map<long int, int> count_map(std::vector<long int> input)
{
    // Define a type for convenience
    typedef std::map<long int,int> CounterMap;
    CounterMap counts;

    // Iterate over our barcodes and increase/initiate counts
    for (size_t i = 0; i < input.size(); i++)
    {
        CounterMap::iterator j(counts.find(input[i]));
        if (j != counts.end())
        {
            j->second++;
        }
        else
        {
            counts[input[i]] = 1;
        }
    }
    return counts;
}

// Helper method - return a count map of mutations in a cell population
std::map<int, int> count_map(std::vector<int> input)
{
    // Define a type for convenience
    typedef std::map<int,int> CounterMap;
    CounterMap counts;

    // Iterate over our barcodes and increase/initiate counts
    for (size_t i = 0; i < input.size(); i++)
    {
        CounterMap::iterator j(counts.find(input[i]));
        if (j != counts.end())
        {
            j->second++;
        }
        else
        {
            counts[input[i]] = 1;
        }
    }
    return counts;
}

// Helper Method - prints a vector of long ints
void printVec(std::vector<long int> v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<long int>(std::cout, " "));
}

// Helper Method - prints a vector of doubles
void printVec(std::vector<double> v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<long int>(std::cout, " "));
}

// Helper method - find the beginning index for a certain triangle and level
int beginIndex(int mfac, int level)
{
    if(level < 2)
    {
        return 0;
    }
    int ret = 0;
    for(int i = 2; i <= level; i++)
    {
        ret = ret + mfac * std::pow(2, i - 2);
    }
    return ret;
}

// Helper method - doubles and splits a vector into groups according to probabilities
std::vector<std::vector<long int>> splitDoubleVector(std::vector<long int> v, std::vector<double> probs)
{
    std::vector<std::vector<long int>> ret = std::vector<std::vector<long int>>(probs.size());

    for(size_t i = 0; i < v.size(); i++)
    {
        int choice = choose(probs);
        ret[choice].push_back(v[i]);
        ret[choice].push_back(v[i]);
    }
    return ret;
}




