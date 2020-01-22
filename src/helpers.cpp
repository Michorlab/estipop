/*
 * =====================================================================================
 *
 *       Filename:  helpers.cpp
 *
 *    Description:  implementing various helper methods
 *
 *        Version:  1.0
 *        Created:  06/13/2018 14:57:27
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
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

#include <RcppGSL.h>
#include <Rcpp.h>
#include <Rinternals.h>

extern gsl_rng* rng;


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

// Maximize a function
double maximizeFunc(gsl_function rate_function, double start_time, double end_time, int bins)
{
  double max = GSL_FN_EVAL(&(rate_function), start_time);
  double test_max = 0;
  double delta_t = (end_time - start_time) / bins;

  if(delta_t > 0.1)
  {
    delta_t = 0.1;
    bins = ceil((end_time - start_time) / delta_t);
  }

  for(int step = 1; step <= bins; ++step)
  {
    test_max = GSL_FN_EVAL(&(rate_function), start_time + delta_t * step);
    max = fmax(max, test_max);
  }
  return max;
}

// Generate piecewise max function for simulation of inhomogenous processes
std::vector<double> maximizePiecewise(gsl_function rate_function, double start_time, double end_time, int bins)
{
  double max = GSL_FN_EVAL(&(rate_function), end_time);
  double test_max = 0;
  double delta_t = (end_time - start_time) / bins;

  std::vector<double> maxes = std::vector<double>(bins);

  for(int step = bins-1; step >= 0; --step)
  {
    test_max = GSL_FN_EVAL(&(rate_function), start_time + delta_t * step);
    max = fmax(max, test_max);
    maxes[step] = max; 
  }
  return maxes;
}

