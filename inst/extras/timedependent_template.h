#ifndef CDIST_H
#define CDIST_H

#include <gsl/gsl_randist.h>
#include <math.h>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif
double rate(double time, void* p);
#ifdef __cplusplus
}
#endif
#endif
