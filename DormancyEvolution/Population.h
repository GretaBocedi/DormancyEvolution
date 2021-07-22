#pragma once

#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>

#include "Individual.h"

using namespace std;

class Population {
public:

	Population();
	~Population();
	int N; //number of adults
	int Nseeds; //number of newly produced seeds
	int NseedsB; //number of seeds in the seedbank 
	double mean_dorm, std_dorm;
	vector<Individual> inds; //adults
	vector<Individual> seeds; //newly produced seeds
	vector<Individual> seedsB; //seedbank
	vector<Individual> seedsB_tmp;

	void reproduce(int, int, int, bool, std::bernoulli_distribution, std::uniform_real_distribution<>);
	void outPop(int, int, int, int, std::ofstream*);
	void outTrait(int, int, int, int, std::ofstream*);

private:
};

