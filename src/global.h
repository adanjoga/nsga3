#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <limits>
#include <algorithm>

#define MAX_INT 	numeric_limits<size_t>::max()
#define MAX_DOUBLE 	numeric_limits<double>::max()
#define MIN_DOUBLE 	numeric_limits<double>::min()

using namespace std;

char strTestInstance[256];

// Bounds of variables
double lowBound = 0, uppBound = 1;

// Dimensionality of variables and objectives
int numVariables;
int numObjectives;

// maximal number of generations
int max_gen;

// Distribution indexes in SBX and polynomial mutation
int id_cx = 30;    // crossover
int id_mu = 20;    // for mutation

// Ideal point used for normalization
double *idealpoint;

// Nadir point used normalization
double *nadirpoint;

// number of reference points for two-objective problems
size_t nrefpoints_2D = 99;

// Number of divisions points in the boundary and inside layers
size_t p_boundary = 0, p_inside = 0;

// Parameters for random number generation
int seed = 237;
long rnd_uni_init;

#include "random.h"

#endif /* GLOBAL_H_ */
