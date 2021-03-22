#pragma once

#include <random>

class Individual {
public:
	Individual(double, int, int);
	~Individual();
	bool alive;
	int x, y;
	double dorm;

	void mutate(std::uniform_real_distribution<>);
	void disp(int, int);
	bool germinate(void);

private:
};

