#pragma once


#include <stdio.h>
#include <stdlib.h>
#include <tchar.h> 
#include <io.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <random>
#include <iterator>

using namespace std;

const int x_max = 200;
const int y_max = 50;

bool fixed_s = true;

int SimNr = 2;
int experiment = 12;
int replicates = 10;
int gen = 500; //generations
int genB = 10;
int out_interval = 5;
int init_x = 2;
int init_y = 25;
int K = 5; //carrying capacity
int s = 5; //nr. seeds per plant
int gamma = 10; //max. nr. of seeds per plant
double d = 0.1; //dispersal probaibility
double m = 0.05; //annual mortality rate of dormant seeds
double beta = 0.01; //mutation probability
double init_dorm = 0.1; //initial dormancy value
double min_fix_dorm = 0.0; //mininimum fix dormancy
double max_fix_dorm = 0.95; //maximum fix dormancy

int Ntot, NseedsB_tot;

ofstream pop, ind;

//-----------------------------------------------------------------------------
// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());

std::uniform_int_distribution<> step(-1, 1);
std::uniform_int_distribution<> sample_s(0, gamma);
std::uniform_real_distribution<> mutation(-0.1, 0.1);
std::bernoulli_distribution disperse(d);
std::bernoulli_distribution dorm_surv(1.0 - m);
std::bernoulli_distribution mutate(beta);



class Individual {
public:
	Individual(double, int, int);
	bool alive;
	int x, y;
	double dorm;

	void mutate(void);
	void disp(int, int);
	bool germinate(void);

private:
};

Individual::Individual(double d, int xx, int yy) {
	alive = true;
	dorm = d;
	x = xx;
	y = yy;
}

class Population {
public:

	Population();
	int N;
	int Nseeds, NseedsB;
	double mean_dorm, std_dorm;
	vector<Individual> inds;
	vector<Individual> seeds, seedsB, seedsB_tmp;

	void reproduce(int, int);
	void outPop(int, int, int, int);
	
private:
};

Population::Population() {
	N = 0;
	Nseeds = 0;
	NseedsB = 0;
	mean_dorm = 0.0;
	std_dorm = 0.0;
}

Population land[x_max][y_max]; 
Population land0[50][50];





//Declaring functions
void Initialise(void);
void InitialiseB(int);
void reproduction(void);
void reproductionB(void);
void dispersal(void);
void dispersalB(void);
void survival(int, int);
void survivalB(int, int);
void outPop_header(double);
void outInd_header(void);
void delete_landscape(void);