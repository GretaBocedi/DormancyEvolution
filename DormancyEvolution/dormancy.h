#pragma once

#define CLUSTER 0

#include <stdio.h>
#include <stdlib.h>
#if CLUSTER 
#include <unistd.h>
#else
#include <tchar.h> 
#include <direct.h>
#include <io.h>
#endif
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <random>
#include <iterator>

#include "Population.h"

using namespace std;

const int x_max = 400;
const int y_max = 50;
const int x_max0 = 50;

bool dorm_evol = false;
bool Allee = false;
bool stochasticity = false;
bool limit_outputs = false;

int SimNr = 2;
int replicates = 10;
int years = 1000; //generations
int expansion_start = 0; //before expansions start dispersal is restricted in 50 x 50 grid
int out_interval = 1;
int outT_interval = 5;
int init_x = 2;
int init_y = 25;
double K = 25.0; //carrying capacity
int s = 5; //nr. seeds per plant
int gam = 10; //max. nr. of seeds per plant
double d = 0.1; //dispersal probaibility
double m = 0.05; //annual mortality rate of dormant seeds
double beta = 0.01; //mutation probability
double init_dorm = 0.1; //initial dormancy value
double env_std = 0.2; //environmental standard deviation
double env_ac = 0.0; //temporal autocorrelation

int Ntot, NseedsB_tot;
int x_m; //changing max x coordinate
double fec = s;
double eps = 0.0;

string dir, dirOut;
ofstream pop, para, trait;

Population land[x_max][y_max];



//-----------------------------------------------------------------------------
// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());



//Declaring functions
void RunModel(void);
void Initialise(void);
void env_stoch(std::normal_distribution<> normEnv);
void reproduction(std::bernoulli_distribution, std::uniform_real_distribution<>);
void dispersal(std::bernoulli_distribution);
void survival(int, int, std::bernoulli_distribution);
void outPara(void);
void outPop_header(void);
void outTrait_header(void);
void delete_landscape(void);